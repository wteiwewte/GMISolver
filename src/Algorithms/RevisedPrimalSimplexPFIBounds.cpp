#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"

#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/SpdlogHeader.h"

namespace {
template <typename T>
void removeElements(std::vector<T> &vec,
                    const std::vector<bool> &shouldBeRemoved) {
  vec.erase(std::remove_if(std::begin(vec), std::end(vec),
                           [&](const T &elem) {
                             return shouldBeRemoved[&elem - &vec[0]];
                           }),
            std::end(vec));
}
} // namespace

template <typename T, typename ComparisonTraitsT>
RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::
    RevisedPrimalSimplexPFIBounds(
        SimplexTableau<T> &simplexTableau,
        const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
        const int32_t objValueLoggingFrequency,
        const int32_t reinversionFrequency)
    : _simplexTableau(simplexTableau),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency) {
  _simplexTableau.calculateRHS();
  _simplexTableau.calculateCurrentObjectiveValue();
  _simplexTableau.calculateSolution();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::run() {
  if (!runPhaseOne()) {
    SPDLOG_WARN("PHASE ONE OF PRIMAL SIMPLEX ALGORITHM FAILED");
    return;
  }

  runPhaseTwo();
}

template <typename T, typename ComparisonTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runPhaseOne() {
  SPDLOG_INFO("BASIS SIZE {}, COLUMN PIVOT RULE {}",
              _simplexTableau._rowInfos.size(),
              primalSimplexColumnPivotRuleToStr(_primalSimplexColumnPivotRule));
  runImpl();

  if (ComparisonTraitsT::greater(_simplexTableau._objectiveValue, 0.0)) {
    SPDLOG_WARN("PROGRAM WITH ARTIFICIAL VARIABLE HAS OPTIMUM GREATER THAN 0 - "
                "INITIAL PROGRAM IS INFEASIBLE");
    return false;
  }

  if (!removeArtificialVariablesFromBasis()) {
    SPDLOG_WARN("COULD NOT REMOVE PROPERLY ARTIFICIAL VARIABLES FROM BASIS");
    return false;
  }
  removeArtificialVariablesFromProgram();
  return true;
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runPhaseTwo() {
  setInitialObjective();
  runImpl();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runImpl() {
  [[maybe_unused]] int iterCount = 1;
  SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  while (true) {
    const bool iterResult = runOneIteration();
    if (iterResult)
      break;

    _simplexTableau.calculateCurrentObjectiveValue();
    _simplexTableau.calculateSolution();

    ++iterCount;
    if (_objValueLoggingFrequency &&
        (iterCount % _objValueLoggingFrequency == 0)) {
      SPDLOG_INFO("ITERATION {}", iterCount);
      SPDLOG_INFO("{}\n", _simplexTableau.toStringObjectiveValue());
    }
    if (_reinversionFrequency && (iterCount % _reinversionFrequency == 0)) {
      if (!reinversion()) {
        SPDLOG_WARN("STOPPING PRIMAL SIMPLEX BECAUSE OF FAILED REINVERSION");
        _simplexTableau._result = LPOptimizationResult::FAILED_REINVERSION;
        break;
      }
    }

    SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  }
  SPDLOG_INFO("PRIMAL SIMPLEX ENDED, ITERATION COUNT {}", iterCount);
}
template <typename T, typename ComparisonTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::runOneIteration() {
  const std::optional<int> enteringColumnIdx = chooseEnteringColumn();
  if (!enteringColumnIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

  SPDLOG_DEBUG("ENTERING COLUMN IDX {} REDUCED COST {}", *enteringColumnIdx,
               _simplexTableau._reducedCosts[*enteringColumnIdx]);

  const std::vector<T> enteringColumn =
      _simplexTableau.computeTableauColumn(*enteringColumnIdx);
  const std::optional<PivotRowData<T>> pivotRowData =
      chooseRowIdx(*enteringColumnIdx, enteringColumn);
  if (!pivotRowData.has_value()) {
    _simplexTableau._result = LPOptimizationResult::UNBOUNDED;
    return true;
  }

  changeTableau(*pivotRowData, *enteringColumnIdx, enteringColumn);
  return false;
}

template <typename T, typename ComparisonTraitsT>
std::optional<int>
RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::chooseEnteringColumn() {
  switch (_primalSimplexColumnPivotRule) {
  case PrimalSimplexColumnPivotRule::FIRST_ELIGIBLE:
    return chooseEnteringColumnFirstEligible();
  case PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST:
    return chooseEnteringColumnBiggestAbsReducedCost();
  }
}

template <typename T, typename ComparisonTraitsT>
std::optional<int> RevisedPrimalSimplexPFIBounds<
    T, ComparisonTraitsT>::chooseEnteringColumnFirstEligible() {
  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size();
       ++columnIdx) {
    SPDLOG_TRACE(
        "COL IDX {} ART {} BASIC {} RED COST {} < 0.0 {} > 0.0 {}", columnIdx,
        _simplexTableau._variableInfos[columnIdx]._isArtificial,
        _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx],
        _simplexTableau._reducedCosts[columnIdx],
        ComparisonTraitsT::less(_simplexTableau._reducedCosts[columnIdx], 0.0),
        ComparisonTraitsT::greater(_simplexTableau._reducedCosts[columnIdx],
                                   0.0));

    SPDLOG_TRACE("COL IDX {}, IS LB {}, IS UB {}", columnIdx,
                 _simplexTableau._simplexBasisData
                     ._isColumnAtLowerBoundBitset[columnIdx],
                 _simplexTableau._simplexBasisData
                     ._isColumnAtUpperBoundBitset[columnIdx]);

    if (!_simplexTableau.isColumnAllowedToEnterBasis(columnIdx))
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
    T, ComparisonTraitsT>::chooseEnteringColumnBiggestAbsReducedCost() {
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
       ++columnIdx)
    if (_simplexTableau.isColumnAllowedToEnterBasis(columnIdx))
      tryUpdateBest(columnIdx);

  return bestColumnIdx;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::changeTableau(
    const PivotRowData<T> &pivotRowData, const int enteringColumnIdx,
    const std::vector<T> &enteringColumn) {
  if (pivotRowData._noBasisChangeNeeded)
    moveVarFromOneBoundToAnother(pivotRowData, enteringColumnIdx,
                                 enteringColumn);
  else
    _simplexTableau.pivotImplicitBounds(
        *pivotRowData._pivotRowIdx, enteringColumnIdx, enteringColumn,
        _simplexTableau.computeTableauRow(*pivotRowData._pivotRowIdx),
        pivotRowData._departingIdxBecomesLowerBound);
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::
    moveVarFromOneBoundToAnother(const PivotRowData<T> &pivotRowData,
                                 const int enteringColumnIdx,
                                 const std::vector<T> &enteringColumn) {
  SPDLOG_DEBUG("NO BASIS CHANGE NEEDED, ENTERING COLUMN IDX {}",
               enteringColumnIdx);

  auto &isColumnAtLowerBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset;
  auto &isColumnAtUpperBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtUpperBoundBitset;
  const bool isVarIncreasing = isColumnAtLowerBoundBitset[enteringColumnIdx];

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    auto addedValue = (enteringColumn[rowIdx] * (*pivotRowData._minRatio));
    if (isVarIncreasing)
      addedValue = -addedValue;

    _simplexTableau._rightHandSides[rowIdx] += addedValue;
  }

  isColumnAtLowerBoundBitset[enteringColumnIdx] =
      !isColumnAtLowerBoundBitset[enteringColumnIdx];
  isColumnAtUpperBoundBitset[enteringColumnIdx] =
      !isColumnAtUpperBoundBitset[enteringColumnIdx];
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
    if (!ComparisonTraitsT::isEligibleForPivot(enteringColumn[rowIdx]))
      continue;

    const auto basicColumnIdx = _simplexTableau.basicColumnIdx(rowIdx);
    const auto &lowerBound =
        _simplexTableau._variableLowerBounds[basicColumnIdx];
    const auto &upperBound =
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
bool RevisedPrimalSimplexPFIBounds<
    T, ComparisonTraitsT>::removeArtificialVariablesFromBasis() {
  std::vector<bool> shouldRowBeRemoved(_simplexTableau._rowInfos.size(), false);

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    SPDLOG_WARN("FOUND BASIC ARTIFICIAL ROW IDX {} COLUMN IDX {}, RHS {}",
                rowIdx, basicVarIdx, _simplexTableau._rightHandSides[rowIdx]);

    std::optional<int> nonZeroEntryColumnIndex;
    const std::vector<T> pivotRow = _simplexTableau.computeTableauRow(rowIdx);

    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
      if (!_simplexTableau.isColumnAllowedToEnterBasis(j))
        continue;

      if (!ComparisonTraitsT::equal(pivotRow[j], 0.0)) {
        nonZeroEntryColumnIndex = j;
        break;
      }
    }

    if (!nonZeroEntryColumnIndex.has_value()) {
      shouldRowBeRemoved[rowIdx] = true;
    } else {
      _simplexTableau.pivotImplicitBounds(
          rowIdx, *nonZeroEntryColumnIndex,
          _simplexTableau.computeTableauColumn(*nonZeroEntryColumnIndex),
          pivotRow,
          true);
    }
  }

  if (std::any_of(shouldRowBeRemoved.begin(), shouldRowBeRemoved.end(),
                  [](const bool val) { return val; })) {
    SPDLOG_INFO("REDUNDANT CONSTRAINTS IN LP FORMULATION");
    removeRows(shouldRowBeRemoved);
    return reinversion();
  }

  return true;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::removeRows(
    const std::vector<bool> &shouldRowBeRemoved) {
  auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset, _1, _2] =
      _simplexTableau._simplexBasisData;

  std::vector<int> oldRowIdxToNewRowIdx(shouldRowBeRemoved.size());
  int curNewRowIdx = 0;

  for (int rowIdx = 0; rowIdx < shouldRowBeRemoved.size(); ++rowIdx)
    if (shouldRowBeRemoved[rowIdx])
      isBasicColumnIndexBitset[rowToBasisColumnIdxMap[rowIdx]] = false;
    else
      oldRowIdxToNewRowIdx[rowIdx] = curNewRowIdx++;

  std::vector<int> newRowToBasisColumnIdxMap(curNewRowIdx);
  for (int rowIdx = 0; rowIdx < shouldRowBeRemoved.size(); ++rowIdx)
    if (!shouldRowBeRemoved[rowIdx])
      newRowToBasisColumnIdxMap[oldRowIdxToNewRowIdx[rowIdx]] =
          rowToBasisColumnIdxMap[rowIdx];

  rowToBasisColumnIdxMap = std::move(newRowToBasisColumnIdxMap);

  removeElements(_simplexTableau._constraintMatrix, shouldRowBeRemoved);
  removeElements(_simplexTableau._rowInfos, shouldRowBeRemoved);
  removeElements(_simplexTableau._rightHandSides, shouldRowBeRemoved);
  removeElements(_simplexTableau._initialRightHandSides, shouldRowBeRemoved);
  for (int i = 0; i < _simplexTableau._basisMatrixInverse.size(); ++i)
    removeElements(_simplexTableau._basisMatrixInverse[i], shouldRowBeRemoved);

  removeElements(_simplexTableau._basisMatrixInverse, shouldRowBeRemoved);
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFIBounds<T,
                                   ComparisonTraitsT>::setInitialObjective() {
  _simplexTableau._objectiveRow =
      _simplexTableau._initialProgram.getObjective();
  _simplexTableau._objectiveRow.resize(_simplexTableau._variableInfos.size());
  calculateDual();
  _simplexTableau.calculateReducedCostsBasedOnDual();
  _simplexTableau.calculateCurrentObjectiveValue();
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
bool RevisedPrimalSimplexPFIBounds<T, ComparisonTraitsT>::reinversion() {
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
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    newBasisMatrixInverse[rowIdx].resize(basisSize);
    newBasisMatrixInverse[rowIdx][rowIdx] = 1.0;
  }

  std::vector<bool> isUnusedColumn(basisSize, true);

  [[maybe_unused]] const auto findPivotColumnFirstEligible =
      [&](const int rowIdx) -> std::optional<int> {
    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx] &&
          !ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
        return colIdx;

    return std::nullopt;
  };

  [[maybe_unused]] const auto findPivotColumnMaxAbsValue =
      [&](const int rowIdx) -> std::optional<int> {
    std::optional<T> maxAbsPivotValue;
    std::optional<int> maxAbsPivotColIdx;
    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx] &&
          ComparisonTraitsT::isEligibleForPivot(basisColumns[colIdx][rowIdx])) {
        const auto currentAbsPivotValue =
            std::fabs(basisColumns[colIdx][rowIdx]);
        if (!maxAbsPivotColIdx.has_value() ||
            (*maxAbsPivotColIdx < currentAbsPivotValue)) {
          maxAbsPivotValue = currentAbsPivotValue;
          maxAbsPivotColIdx = colIdx;
        }
      }

    return maxAbsPivotColIdx;
  };

  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    const auto pivotColumnIdx = findPivotColumnMaxAbsValue(rowIdx);
    if (!pivotColumnIdx.has_value()) {
      SPDLOG_WARN("Basis matrix reinversion failed for row {}!", rowIdx);
      return false;
    }

    const T pivotingTermInverse{1.0 / basisColumns[*pivotColumnIdx][rowIdx]};
    for (int j = 0; j < basisSize; ++j) {
      if (j == rowIdx)
        continue;

      //      if (ComparisonTraitsT::equal(basisColumns[*pivotColumnIdx][j],
      //      0.0))
      //        continue;

      const auto commonCoeff =
          pivotingTermInverse * basisColumns[*pivotColumnIdx][j];

      //      if (ComparisonTraitsT::equal(commonCoeff, 0.0))
      //        continue;
      for (int k = 0; k < basisSize; ++k)
        newBasisMatrixInverse[j][k] -=
            commonCoeff * newBasisMatrixInverse[rowIdx][k];
    }

    for (int k = 0; k < basisSize; ++k)
      newBasisMatrixInverse[rowIdx][k] *= pivotingTermInverse;

    isUnusedColumn[*pivotColumnIdx] = false;
    _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
        columnIndexesMapping[*pivotColumnIdx];

    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx]) {
        //        if (ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx],
        //        0.0))
        //          continue;

        const auto commonCoeff =
            pivotingTermInverse * basisColumns[colIdx][rowIdx];

        //        if (!ComparisonTraitsT::equal(commonCoeff, 0.0))
        for (int k = 0; k < basisSize; ++k) {
          if (k == rowIdx)
            continue;

          basisColumns[colIdx][k] -=
              commonCoeff * basisColumns[*pivotColumnIdx][k];
        }

        basisColumns[colIdx][rowIdx] *= pivotingTermInverse;
      }
  }

  _simplexTableau._basisMatrixInverse.swap(newBasisMatrixInverse);
  _simplexTableau.calculateRHS();
  SPDLOG_INFO("REINVERSION SUCCESS");
  return true;
}

template class RevisedPrimalSimplexPFIBounds<double>;