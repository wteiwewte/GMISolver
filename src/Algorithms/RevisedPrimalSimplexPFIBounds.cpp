#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"

#include "src/Algorithms/SimplexTableau.h"

template <typename T, typename ComparisonTraitsT>
RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::
    RevisedPrimalSimplexPFIBounds(SimplexTableau<T> &simplexTableau)
    : _simplexTableau(simplexTableau) {
  initRHS();
  calculateCurrentObjectiveValue();
  calculateSolution();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::run() {
  if (!runPhaseOne())
    return;

  runPhaseTwo();
}

template <typename T, typename ComparisonTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runPhaseOne() {
  SPDLOG_INFO("BASIS SIZE {}", _simplexTableau._rowInfos.size());
  runImpl();

  if (ComparisonTraitsT::greater(_simplexTableau._objectiveValue, 0.0)) {
    SPDLOG_INFO("Program with artificial variable has optimum greater than 0 - "
                "initial program is infeasible");
    return false;
  }

  return true;
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runPhaseTwo() {
  removeArtificialVariablesFromBasis();
  removeArtificialVariablesFromProgram();
  setInitialObjective();
  runImpl();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runImpl() {
  [[maybe_unused]] int iterCount = 1;
  SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  while (true) {
    SPDLOG_DEBUG("ITERATION {}", iterCount++);
    const bool iterResult = runOneIteration();
    if (iterResult)
      break;

    calculateCurrentObjectiveValue();
    calculateSolution();

    SPDLOG_DEBUG("{}\n", _simplexTableau.toStringShort());
    SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  }
  SPDLOG_INFO("PRIMAL SIMPLEX ENDED");
}
template <typename T, typename ComparisonTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runOneIteration() {
  const std::optional<int> enteringColumnIdx =
      chooseEnteringColumnIdxBiggestReducedCost();
  if (!enteringColumnIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

  SPDLOG_DEBUG("ENTERING COLUMN IDX {} REDUCED COST {}", *enteringColumnIdx,
               _simplexTableau._reducedCosts[*enteringColumnIdx]);

  const std::vector<T> enteringColumn =
      computeEnteringColumn(*enteringColumnIdx);
  const std::optional<PivotRowData<T>> pivotRowData =
      chooseRowIdx(*enteringColumnIdx, enteringColumn);
  if (!pivotRowData.has_value()) {
    _simplexTableau._result = LPOptimizationResult::UNBOUNDED;
    return true;
  }

  pivot(*pivotRowData, *enteringColumnIdx, enteringColumn);
  return false;
}

template <typename T, typename ComparisonTraitsT>
std::optional<int>
RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::chooseEnteringColumnIdx() {
  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size();
       ++columnIdx) {
    SPDLOG_DEBUG(
        "COL IDX {} ART {} BASIC {} RED COST {} < 0.0 {} > 0.0 {}", columnIdx,
        _simplexTableau._variableInfos[columnIdx]._isArtificial,
        _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx],
        _simplexTableau._reducedCosts[columnIdx],
        ComparisonTraitsT::less(_simplexTableau._reducedCosts[columnIdx], 0.0),
        ComparisonTraitsT::greater(_simplexTableau._reducedCosts[columnIdx],
                                   0.0));

    SPDLOG_DEBUG("COL IDX {}, IS LB {}, IS UB {}", columnIdx,
                 _simplexTableau._simplexBasisData
                     ._isColumnAtLowerBoundBitset[columnIdx],
                 _simplexTableau._simplexBasisData
                     ._isColumnAtUpperBoundBitset[columnIdx]);

    if (_simplexTableau._variableInfos[columnIdx]._isArtificial)
      continue;

    if (_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx])
      continue;

    if (ComparisonTraitsT::less(_simplexTableau._reducedCosts[columnIdx],
                                0.0) &&
        _simplexTableau._simplexBasisData
            ._isColumnAtLowerBoundBitset[columnIdx])
      return columnIdx;

    if (ComparisonTraitsT::greater(_simplexTableau._reducedCosts[columnIdx],
                                   0.0) &&
        _simplexTableau._simplexBasisData
            ._isColumnAtUpperBoundBitset[columnIdx])
      return columnIdx;
  }
  return std::nullopt;
}

template <typename T, typename ComparisonTraitsT>
std::optional<int> RevisedPrimalSimplexPFIBounds<
    T, ComparisonTraitsT>::chooseEnteringColumnIdxBiggestReducedCost() {
  std::optional<int> bestColumnIdx;
  std::optional<T> biggestAbsReducedCost;

  const auto tryUpdateBest = [&](const int colIdx) {
    if ((_simplexTableau._simplexBasisData
             ._isColumnAtLowerBoundBitset[colIdx] &&
         ComparisonTraitsT::less(_simplexTableau._reducedCosts[colIdx], 0.0)) ||
        (_simplexTableau._simplexBasisData
             ._isColumnAtUpperBoundBitset[colIdx] &&
         ComparisonTraitsT::greater(_simplexTableau._reducedCosts[colIdx],
                                    0.0))) {
      const auto currentAbsReducedCost =
          std::fabs(_simplexTableau._reducedCosts[colIdx]);
      if (!biggestAbsReducedCost.has_value() ||
          (*biggestAbsReducedCost < currentAbsReducedCost)) {
        bestColumnIdx = colIdx;
        biggestAbsReducedCost = currentAbsReducedCost;
      }
    }
  };

  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size();
       ++columnIdx) {
    if (_simplexTableau._variableInfos[columnIdx]._isArtificial)
      continue;

    if (_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx])
      continue;

    tryUpdateBest(columnIdx);
  }
  return bestColumnIdx;
}

template <typename T, typename ComparisonTraitsT>
std::vector<T>
RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::computeEnteringColumn(
    const int enteringColumnIdx) {
  std::vector<T> result(_simplexTableau._rowInfos.size());

  // TODO - maybe add kahan summation algo, maybe opt order
  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i) {
    for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k)
      result[i] += _simplexTableau._basisMatrixInverse[i][k] *
                   _simplexTableau._constraintMatrix[k][enteringColumnIdx];
  }

  return result;
}

template <typename T, typename ComparisonTraitsT>
std::vector<T>
RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::computePivotRow(
    const int rowIdx) {
  std::vector<T> result(_simplexTableau._variableInfos.size());

  // TODO - maybe add kahan summation algo, maybe opt order
  const int columnBasicIdx =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
    if (j == columnBasicIdx)
      result[j] = 1.0;
    else if (j != columnBasicIdx &&
             _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
      result[j] = 0.0;
    else
      for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k)
        result[j] += _simplexTableau._basisMatrixInverse[rowIdx][k] *
                     _simplexTableau._constraintMatrix[k][j];
  }

  return result;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::pivot(
    const PivotRowData<T> &pivotRowData, const int enteringColumnIdx,
    const std::vector<T> &enteringColumn) {
  auto &isColumnAtLowerBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset;
  auto &isColumnAtUpperBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtUpperBoundBitset;
  const bool isVarIncreasing = isColumnAtLowerBoundBitset[enteringColumnIdx];

  if (pivotRowData._noBasisChangeNeeded) {
    SPDLOG_DEBUG("NO BASIS CHANGE NEEDED, ENTERING COLUMN IDX {}",
                 enteringColumnIdx);
    //    isColumnAtLowerBoundBitset[enteringColumnIdx] =
    //    !isColumnAtLowerBoundBitset[enteringColumnIdx];
    //    isColumnAtUpperBoundBitset[enteringColumnIdx] =
    //    !isColumnAtUpperBoundBitset[enteringColumnIdx];

    for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
      auto addedValue = (enteringColumn[rowIdx] * (*pivotRowData._minRatio));
      if (isVarIncreasing)
        addedValue = -addedValue;

      _simplexTableau._rightHandSides[rowIdx] += addedValue;
    }

    // reduced costs?

    isColumnAtLowerBoundBitset[enteringColumnIdx] =
        !isColumnAtLowerBoundBitset[enteringColumnIdx];
    isColumnAtUpperBoundBitset[enteringColumnIdx] =
        !isColumnAtUpperBoundBitset[enteringColumnIdx];
  } else {
    SPDLOG_DEBUG("PIVOT - ENTERING COLUMN IDX {}, ROW IDX {}",
                 enteringColumnIdx, *pivotRowData._pivotRowIdx);
    const auto leavingBasicColumnIdx =
        _simplexTableau._simplexBasisData
            ._rowToBasisColumnIdxMap[*pivotRowData._pivotRowIdx];
    const auto leavingColumn = computeEnteringColumn(leavingBasicColumnIdx);

    for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
      _simplexTableau._rightHandSides[rowIdx] +=
          (isColumnAtLowerBoundBitset[enteringColumnIdx]
               ? *_simplexTableau._variableLowerBounds[enteringColumnIdx]
               : *_simplexTableau._variableUpperBounds[enteringColumnIdx]) *
          enteringColumn[rowIdx];

      _simplexTableau._rightHandSides[rowIdx] -=
          (pivotRowData._departingIdxBecomesLowerBound
               ? *_simplexTableau._variableLowerBounds[leavingBasicColumnIdx]
               : *_simplexTableau._variableUpperBounds[leavingBasicColumnIdx]) *
          leavingColumn[rowIdx];
    }

    executePivot(*pivotRowData._pivotRowIdx, enteringColumnIdx, enteringColumn,
                 computePivotRow(*pivotRowData._pivotRowIdx));
    (pivotRowData._departingIdxBecomesLowerBound
         ? isColumnAtLowerBoundBitset[leavingBasicColumnIdx]
         : isColumnAtUpperBoundBitset[leavingBasicColumnIdx]) = true;

    if (isColumnAtLowerBoundBitset[enteringColumnIdx])
      isColumnAtLowerBoundBitset[enteringColumnIdx] = false;

    if (isColumnAtUpperBoundBitset[enteringColumnIdx])
      isColumnAtUpperBoundBitset[enteringColumnIdx] = false;
  }
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::executePivot(
    const int rowIdx, const int enteringColumnIdx,
    const std::vector<T> &enteringColumn, const std::vector<T> &pivotRow) {
  const PivotData<T> pivotData{rowIdx, enteringColumnIdx,
                               1.0 / enteringColumn[rowIdx]};
  updateReducedCosts(pivotData, computePivotRow(rowIdx));
  SPDLOG_DEBUG("PIVOT VALUE {},  INV {}", enteringColumn[rowIdx],
               1.0 / enteringColumn[rowIdx]);
  updateInverseMatrixWithRHS(pivotData, enteringColumn);
  //    calculateDual();
  //    _simplexTableau.calculateReducedCostsBasedOnDual();
  _simplexTableau.updateBasisData(pivotData);
}

template <typename T, typename ComparisonTraitsT>
std::optional<PivotRowData<T>>
RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::chooseRowIdx(
    const int enteringColumnIdx, const std::vector<T> &enteringColumn) {
  std::optional<int> rowIdxMostRestrictiveLowerBound;
  std::optional<T> mostRestrictiveLowerBound;

  std::optional<int> rowIdxMostRestrictiveUpperBound;
  std::optional<T> mostRestrictiveUpperBound;

  for (int rowIdx = 0; rowIdx < _simplexTableau.getRowInfos().size();
       ++rowIdx) {
    if (ComparisonTraitsT::equal(enteringColumn[rowIdx], 0.0))
      continue;

    const auto basicColumnIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    const auto lowerBound =
        _simplexTableau._variableLowerBounds[basicColumnIdx];
    const auto upperBound =
        _simplexTableau._variableUpperBounds[basicColumnIdx];

    if (_simplexTableau._simplexBasisData
            ._isColumnAtLowerBoundBitset[enteringColumnIdx]) {
      if (enteringColumn[rowIdx] > 0.0) {
        if (!lowerBound.has_value())
          continue;

        const auto currentRestrictionValue =
            (_simplexTableau._rightHandSides[rowIdx] - *lowerBound) /
            enteringColumn[rowIdx];
        if (!rowIdxMostRestrictiveLowerBound.has_value() ||
            (currentRestrictionValue < *mostRestrictiveLowerBound)) {
          rowIdxMostRestrictiveLowerBound = rowIdx;
          mostRestrictiveLowerBound = currentRestrictionValue;
        }
      } else {
        if (!upperBound.has_value())
          continue;

        const auto currentRestrictionValue =
            (*upperBound - _simplexTableau._rightHandSides[rowIdx]) /
            (-enteringColumn[rowIdx]);
        if (!rowIdxMostRestrictiveUpperBound.has_value() ||
            (currentRestrictionValue < *mostRestrictiveUpperBound)) {
          rowIdxMostRestrictiveUpperBound = rowIdx;
          mostRestrictiveUpperBound = currentRestrictionValue;
        }
      }
    } else {
      if (enteringColumn[rowIdx] > 0.0) {
        if (!upperBound.has_value())
          continue;

        const auto currentRestrictionValue =
            (*upperBound - _simplexTableau._rightHandSides[rowIdx]) /
            enteringColumn[rowIdx];
        if (!rowIdxMostRestrictiveUpperBound.has_value() ||
            (currentRestrictionValue < *mostRestrictiveUpperBound)) {
          rowIdxMostRestrictiveUpperBound = rowIdx;
          mostRestrictiveUpperBound = currentRestrictionValue;
        }
      } else {
        if (!lowerBound.has_value())
          continue;

        const auto currentRestrictionValue =
            (_simplexTableau._rightHandSides[rowIdx] - *lowerBound) /
            (-enteringColumn[rowIdx]);
        if (!rowIdxMostRestrictiveLowerBound.has_value() ||
            (currentRestrictionValue < *mostRestrictiveLowerBound)) {
          rowIdxMostRestrictiveLowerBound = rowIdx;
          mostRestrictiveLowerBound = currentRestrictionValue;
        }
      }
    }
  }

  PivotRowData<T> pivotRowData;
  const auto enteringColumnLowerBound =
      _simplexTableau._variableLowerBounds[enteringColumnIdx];
  const auto enteringColumnUpperBound =
      _simplexTableau._variableUpperBounds[enteringColumnIdx];
  if (enteringColumnLowerBound.has_value() &&
      enteringColumnUpperBound.has_value()) {
    pivotRowData._minRatio =
        *enteringColumnUpperBound - *enteringColumnLowerBound;
    pivotRowData._noBasisChangeNeeded = true;
  }

  if (rowIdxMostRestrictiveLowerBound.has_value()) {
    if (!pivotRowData._minRatio.has_value() ||
        (*mostRestrictiveLowerBound < pivotRowData._minRatio)) {
      pivotRowData._minRatio = *mostRestrictiveLowerBound;
      pivotRowData._pivotRowIdx = *rowIdxMostRestrictiveLowerBound;
      pivotRowData._departingIdxBecomesLowerBound = true;
      pivotRowData._noBasisChangeNeeded = false;
    }
  }

  if (rowIdxMostRestrictiveUpperBound.has_value()) {
    if (!pivotRowData._minRatio.has_value() ||
        (*mostRestrictiveUpperBound < pivotRowData._minRatio)) {
      pivotRowData._minRatio = *mostRestrictiveUpperBound;
      pivotRowData._pivotRowIdx = *rowIdxMostRestrictiveUpperBound;
      pivotRowData._departingIdxBecomesLowerBound = false;
      pivotRowData._noBasisChangeNeeded = false;
    }
  }

  return pivotRowData._minRatio.has_value() ? std::optional{pivotRowData}
                                            : std::nullopt;
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::updateReducedCosts(
    const PivotData<T> &pivotData, const std::vector<T> &pivotRow) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;
  const auto commonCoeffReducedCost =
      _simplexTableau._reducedCosts[enteringColumnIdx] * pivotingTermInverse;
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    _simplexTableau._reducedCosts[j] -= commonCoeffReducedCost * pivotRow[j];
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::
    updateInverseMatrixWithRHS(const PivotData<T> &pivotData,
                               const std::vector<T> &enteringColumn) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;

  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i) {
    if (i == leavingRowIdx)
      continue;
    // TODO - opt if coeff is zero
    const auto commonCoeff = enteringColumn[i] * pivotingTermInverse;

    for (int j = 0; j < _simplexTableau._basisMatrixInverse.size(); ++j)
      _simplexTableau._basisMatrixInverse[i][j] -=
          commonCoeff * _simplexTableau._basisMatrixInverse[leavingRowIdx][j];

    _simplexTableau._rightHandSides[i] -=
        commonCoeff * _simplexTableau._rightHandSides[leavingRowIdx];
  }

  for (int j = 0; j < _simplexTableau._basisMatrixInverse.size(); ++j)
    _simplexTableau._basisMatrixInverse[leavingRowIdx][j] *=
        pivotingTermInverse;

  _simplexTableau._rightHandSides[leavingRowIdx] *= pivotingTermInverse;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<
    T, ComparisonTraitsT>::removeArtificialVariablesFromProgram() {
  std::optional<int> firstArtificialIdx;
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    if (_simplexTableau._variableInfos[j]._isArtificial) {
      firstArtificialIdx = j;
      break;
    }

  if (!firstArtificialIdx.has_value())
    return;

  _simplexTableau._variableInfos.resize(*firstArtificialIdx);
  _simplexTableau._reducedCosts.resize(*firstArtificialIdx);
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    _simplexTableau._constraintMatrix[rowIdx].resize(*firstArtificialIdx);
  }
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<
    T, ComparisonTraitsT>::removeArtificialVariablesFromBasis() {
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicVarIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    SPDLOG_INFO("FOUND BASIC ARTIFICIAL ROW IDX {} COLUMN IDX {}", rowIdx,
                basicVarIdx);

    std::optional<int> nonZeroEntryColumnIndex;
    const std::vector<T> pivotRow = computePivotRow(rowIdx);

    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
      if (_simplexTableau._variableInfos[j]._isArtificial ||
          _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        continue;

      if (!ComparisonTraitsT::equal(pivotRow[j], 0.0)) {
        nonZeroEntryColumnIndex = j;
        break;
      }
    }

    if (!nonZeroEntryColumnIndex.has_value()) {
      SPDLOG_WARN("Redundant constraints in lp formulation!");
      removeRow(rowIdx);
    } else {
      // FIXME
      SPDLOG_INFO("{} ENTERING COLUMN IDX", *nonZeroEntryColumnIndex);
      const std::vector<T> enteringColumn =
          computeEnteringColumn(*nonZeroEntryColumnIndex);
      executePivot(rowIdx, *nonZeroEntryColumnIndex, enteringColumn, pivotRow);
    }
  }
  reinversion();
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::removeRow(
    const int rowIdx) {
  // rethink it
  // FIXME - case when multiple rows are deleted air05 case
  auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset, _1, _2] =
      _simplexTableau._simplexBasisData;
  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[rowIdx]] = false;
  for (int i = rowIdx; i < _simplexTableau._rowInfos.size() - 1; ++i)
    rowToBasisColumnIdxMap[i] = rowToBasisColumnIdxMap[i + 1];
  rowToBasisColumnIdxMap.pop_back();

  _simplexTableau._constraintMatrix.erase(
      _simplexTableau._constraintMatrix.begin() + rowIdx);
  _simplexTableau._rowInfos.erase(_simplexTableau._rowInfos.begin() + rowIdx);
  _simplexTableau._rightHandSides.erase(
      _simplexTableau._rightHandSides.begin() + rowIdx);
  _simplexTableau._initialRightHandSides.erase(
      _simplexTableau._initialRightHandSides.begin() + rowIdx);
  for (int i = 0; i < _simplexTableau._basisMatrixInverse.size(); ++i)
    _simplexTableau._basisMatrixInverse[i].erase(
        _simplexTableau._basisMatrixInverse[i].begin() + rowIdx);

  _simplexTableau._basisMatrixInverse.erase(
      _simplexTableau._basisMatrixInverse.begin() + rowIdx);
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T,
                                   ComparisonTraitsT>::setInitialObjective() {
  _simplexTableau._objectiveRow =
      _simplexTableau._initialProgram.getObjective();
  _simplexTableau._objectiveRow.resize(_simplexTableau._variableInfos.size());
  calculateDual();
  _simplexTableau.calculateReducedCostsBasedOnDual();
  calculateCurrentObjectiveValue();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::calculateDual() {
  _simplexTableau._y.resize(_simplexTableau._basisMatrixInverse.size());
  for (int colIdx = 0; colIdx < _simplexTableau._y.size(); ++colIdx) {
    T sum{};

    for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k) {
      SPDLOG_TRACE("{} {} {}", k,
                   _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[k],
                   _simplexTableau._objectiveRow.size());
      sum += _simplexTableau._objectiveRow[_simplexTableau._simplexBasisData
                                               ._rowToBasisColumnIdxMap[k]] *
             _simplexTableau._basisMatrixInverse[k][colIdx];
    }

    _simplexTableau._y[colIdx] = sum;
  }
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<
    T, ComparisonTraitsT>::calculateCurrentObjectiveValue() {
  _simplexTableau._objectiveValue = T{};
  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i)
    _simplexTableau._objectiveValue +=
        _simplexTableau._rightHandSides[i] *
        _simplexTableau._objectiveRow[_simplexTableau._simplexBasisData
                                          ._rowToBasisColumnIdxMap[i]];

  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
      _simplexTableau._objectiveValue +=
          (_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[j]
               ? *_simplexTableau._variableLowerBounds[j]
               : *_simplexTableau._variableUpperBounds[j]) *
          _simplexTableau._objectiveRow[j];
}

// template <typename T, typename ComparisonTraitsT>
// void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::calculateSolution()
// {
//   _simplexTableau._x.resize(_simplexTableau._variableInfos.size());
//   std::fill(_simplexTableau._x.begin(), _simplexTableau._x.end(), T{});
//   for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
//     _simplexTableau._x[_simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx]]
//     =
//         _simplexTableau._rightHandSides[rowIdx];
//
//   for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
//     if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
//       _simplexTableau._x[j]  =
//       (_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[j] ?
//                                                                                                 *_simplexTableau._variableLowerBounds[j] : *_simplexTableau._variableUpperBounds[j]);
// }
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::calculateSolution() {
  _simplexTableau._x.resize(_simplexTableau._variableInfos.size());
  std::fill(_simplexTableau._x.begin(), _simplexTableau._x.end(), T{});

  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
      _simplexTableau._x[j] =
          (_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[j]
               ? *_simplexTableau._variableLowerBounds[j]
               : *_simplexTableau._variableUpperBounds[j]);

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicColumnIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    _simplexTableau._x[basicColumnIdx] =
        _simplexTableau._rightHandSides[rowIdx];

    //    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    //      if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
    //        _simplexTableau._x[basicColumnIdx] -= _simplexTableau._x[j] *
    //        _simplexTableau._constraintMatrix[rowIdx][j];
  }
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::initRHS() {
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    //    const auto basicColumnIdx =
    //    _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
      if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        _simplexTableau._rightHandSides[rowIdx] -=
            (_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[j]
                 ? *_simplexTableau._variableLowerBounds[j]
                 : *_simplexTableau._variableUpperBounds[j]) *
            _simplexTableau._constraintMatrix[rowIdx][j];
  }
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::calculateRHS() {
  std::vector<T> tempRHS = _simplexTableau._initialRightHandSides;
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    //    const auto basicColumnIdx = _simplexTableau._simplexB
    //    asisData._rowToBasisColumnIdxMap[rowIdx];
    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
      if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        tempRHS[rowIdx] -=
            (_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[j]
                 ? *_simplexTableau._variableLowerBounds[j]
                 : *_simplexTableau._variableUpperBounds[j]) *
            _simplexTableau._constraintMatrix[rowIdx][j];
  }

  std::fill(_simplexTableau._rightHandSides.begin(),
            _simplexTableau._rightHandSides.end(), T{0.0});
  // TODO - summation
  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i)
    for (int j = 0; j < _simplexTableau._rowInfos.size(); ++j)
      _simplexTableau._rightHandSides[i] +=
          _simplexTableau._basisMatrixInverse[i][j] * tempRHS[j];
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::reinversion() {
  // TODO - opt it by storing columns in vectors
  const auto basisSize = _simplexTableau._rowInfos.size();
  std::vector<int> columnIndexesMapping(basisSize);
  std::vector<std::vector<T>> basisColumns(basisSize);
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    const auto basicColumnIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    columnIndexesMapping[rowIdx] = basicColumnIdx;
    basisColumns[rowIdx].resize(basisSize);
    for (int j = 0; j < basisSize; ++j)
      basisColumns[rowIdx][j] =
          _simplexTableau._constraintMatrix[j][basicColumnIdx];
  }

  std::vector<std::vector<T>> newBasisMatrixInverse(basisSize);
  std::vector<T> newRHS = _simplexTableau._initialProgram.getRightHandSides();
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    newBasisMatrixInverse[rowIdx].resize(basisSize);
    newBasisMatrixInverse[rowIdx][rowIdx] = 1.0;
  }

  std::vector<bool> isUnusedColumn(basisSize, true);

  const auto findPivotColumn = [&](const int rowIdx) -> std::optional<int> {
    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx] &&
          !ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
        return colIdx;

    return std::nullopt;
  };

  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    const auto pivotColumnIdx = findPivotColumn(rowIdx);
    if (!pivotColumnIdx.has_value()) {
      SPDLOG_WARN("Basis matrix reinversion failed for row {}!", rowIdx);
      return;
    }

    const T pivotingTermInverse{1.0 / basisColumns[*pivotColumnIdx][rowIdx]};
    //    SPDLOG_DEBUG("REINVERSION PIVOT {}", pivotingTermInverse);
    for (int j = 0; j < basisSize; ++j) {
      if (j == rowIdx)
        continue; // !!

      if (ComparisonTraitsT::equal(basisColumns[*pivotColumnIdx][j], 0.0))
        continue;

      const auto commonCoeff =
          pivotingTermInverse * basisColumns[*pivotColumnIdx][j];
      for (int k = 0; k < basisSize; ++k)
        newBasisMatrixInverse[j][k] -=
            commonCoeff * newBasisMatrixInverse[rowIdx][k];

      newRHS[j] -= commonCoeff * newRHS[rowIdx];
    }

    for (int k = 0; k < basisSize; ++k)
      newBasisMatrixInverse[rowIdx][k] *= pivotingTermInverse;

    newRHS[rowIdx] *= pivotingTermInverse;

    isUnusedColumn[*pivotColumnIdx] = false;
    _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
        columnIndexesMapping[*pivotColumnIdx];

    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx]) {
        if (ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
          continue;

        for (int k = 0; k < basisSize; ++k) {
          if (k == rowIdx)
            continue;

          basisColumns[colIdx][k] -= pivotingTermInverse *
                                     basisColumns[*pivotColumnIdx][k] *
                                     basisColumns[colIdx][rowIdx];
        }

        basisColumns[colIdx][rowIdx] *= pivotingTermInverse;
      }
  }

  _simplexTableau._basisMatrixInverse = newBasisMatrixInverse;
  //  _simplexTableau._rightHandSides = newRHS;
  calculateRHS();
  SPDLOG_INFO("REINVERSION SUCCESS");
}
// template <typename T, typename ComparisonTraitsT>
// void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::calculateRHS() {
//   for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
//   {
//
//   }
// }

template class RevisedPrimalSimplexPFIBounds<double>;