#include "src/Algorithms/RevisedPrimalSimplexPFIBoundsSparse.h"

#include "src/Algorithms/SimplexTableau.h"
#include "src/DataModel/LinearProgram.h"
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

template <typename T, typename SimplexTraitsT>
RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::
    RevisedPrimalSimplexPFIBoundsSparse(
        SimplexTableau<T, SimplexTraitsT> &simplexTableau,
        const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
        const int32_t objValueLoggingFrequency,
        const int32_t reinversionFrequency)
    : _simplexTableau(simplexTableau),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency) {
  _simplexTableau.calculateRHSPFISparse();
  _simplexTableau.calculateCurrentObjectiveValue();
  _simplexTableau.calculateSolution();
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::run() {
  if (!runPhaseOne()) {
    SPDLOG_WARN("PHASE ONE OF PRIMAL SIMPLEX ALGORITHM SPARSE FAILED");
    return;
  }

  runPhaseTwo();
}

template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::runPhaseOne() {
  SPDLOG_INFO("BASIS SIZE {}, COLUMN PIVOT RULE {}",
              _simplexTableau._rowInfos.size(),
              primalSimplexColumnPivotRuleToStr(_primalSimplexColumnPivotRule));
  runImpl();

  if (SimplexTraitsT::greater(_simplexTableau._objectiveValue, 0.0)) {
    SPDLOG_WARN(
        "PROGRAM WITH ARTIFICIAL VARIABLE HAS OPTIMUM {} GREATER THAN 0 - "
        "INITIAL PROGRAM IS INFEASIBLE",
        _simplexTableau._objectiveValue);
    return false;
  }

  //  if (!removeArtificialVariablesFromBasis()) {
  //    SPDLOG_WARN("COULD NOT REMOVE PROPERLY ARTIFICIAL VARIABLES FROM
  //    BASIS"); return false;
  //  }
  //  removeArtificialVariablesFromProgram();
  return true;
}
template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::runPhaseTwo() {
  _simplexTableau.setObjectiveSparse(
      _simplexTableau._initialProgram.getObjective());
  runImpl();
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::runImpl() {
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
      if (!_simplexTableau.reinversionPFISparse()) {
        SPDLOG_WARN("STOPPING PRIMAL SIMPLEX BECAUSE OF FAILED REINVERSION");
        _simplexTableau._result = LPOptimizationResult::FAILED_REINVERSION;
        break;
      }
    }

    SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  }
  SPDLOG_INFO("PRIMAL SIMPLEX SPARSE ENDED, ITERATION COUNT {}", iterCount);
}
template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::runOneIteration() {
  const std::optional<int> enteringColumnIdx = chooseEnteringColumn();
  if (!enteringColumnIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

  SPDLOG_DEBUG("ENTERING COLUMN IDX {} REDUCED COST {}", *enteringColumnIdx,
               _simplexTableau._reducedCosts[*enteringColumnIdx]);

  const auto enteringColumn =
      _simplexTableau.computeTableauColumnPFISparse(*enteringColumnIdx);
  const std::optional<PivotRowData<T>> pivotRowData =
      chooseRowIdx(*enteringColumnIdx, enteringColumn._normalVec);
  if (!pivotRowData.has_value()) {
    _simplexTableau._result = LPOptimizationResult::UNBOUNDED;
    return true;
  }

  changeTableau(*pivotRowData, *enteringColumnIdx, enteringColumn);
  return false;
}

template <typename T, typename SimplexTraitsT>
std::optional<int>
RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::chooseEnteringColumn() {
  switch (_primalSimplexColumnPivotRule) {
  case PrimalSimplexColumnPivotRule::FIRST_ELIGIBLE:
    return chooseEnteringColumnFirstEligible();
  case PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST:
    return chooseEnteringColumnBiggestAbsReducedCost();
  }
}

template <typename T, typename SimplexTraitsT>
std::optional<int> RevisedPrimalSimplexPFIBoundsSparse<
    T, SimplexTraitsT>::chooseEnteringColumnFirstEligible() {
  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size();
       ++columnIdx) {
    SPDLOG_TRACE(
        "COL IDX {} ART {} BASIC {} RED COST {} < 0.0 {} > 0.0 {}", columnIdx,
        _simplexTableau._variableInfos[columnIdx]._isArtificial,
        _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx],
        _simplexTableau._reducedCosts[columnIdx],
        SimplexTraitsT::less(_simplexTableau._reducedCosts[columnIdx], 0.0),
        SimplexTraitsT::greater(_simplexTableau._reducedCosts[columnIdx], 0.0));

    SPDLOG_TRACE("COL IDX {}, IS LB {}, IS UB {}", columnIdx,
                 _simplexTableau._simplexBasisData
                     ._isColumnAtLowerBoundBitset[columnIdx],
                 _simplexTableau._simplexBasisData
                     ._isColumnAtUpperBoundBitset[columnIdx]);

    if (!_simplexTableau.isColumnAllowedToEnterBasis(columnIdx))
      continue;

    if (SimplexTraitsT::less(_simplexTableau._reducedCosts[columnIdx], 0.0) &&
        _simplexTableau._simplexBasisData
            ._isColumnAtLowerBoundBitset[columnIdx])
      return columnIdx;

    if (SimplexTraitsT::greater(_simplexTableau._reducedCosts[columnIdx],
                                0.0) &&
        _simplexTableau._simplexBasisData
            ._isColumnAtUpperBoundBitset[columnIdx])
      return columnIdx;
  }
  return std::nullopt;
}

template <typename T, typename SimplexTraitsT>
std::optional<int> RevisedPrimalSimplexPFIBoundsSparse<
    T, SimplexTraitsT>::chooseEnteringColumnBiggestAbsReducedCost() {
  std::optional<int> bestColumnIdx;
  std::optional<T> biggestAbsReducedCost;

  const auto tryUpdateBest = [&](const int colIdx) {
    if ((_simplexTableau._simplexBasisData
             ._isColumnAtLowerBoundBitset[colIdx] &&
         SimplexTraitsT::less(_simplexTableau._reducedCosts[colIdx], 0.0)) ||
        (_simplexTableau._simplexBasisData
             ._isColumnAtUpperBoundBitset[colIdx] &&
         SimplexTraitsT::greater(_simplexTableau._reducedCosts[colIdx], 0.0))) {
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

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::changeTableau(
    const PivotRowData<T> &pivotRowData, const int enteringColumnIdx,
    const SparseVector<T> &enteringColumn) {
  if (pivotRowData._noBasisChangeNeeded)
    moveVarFromOneBoundToAnother(pivotRowData, enteringColumnIdx,
                                 enteringColumn._normalVec);
  else
    _simplexTableau.pivotImplicitBoundsSparse(
        *pivotRowData._pivotRowIdx, enteringColumnIdx, enteringColumn,
        _simplexTableau.computeTableauRowPFISparse(*pivotRowData._pivotRowIdx),
        pivotRowData._departingIdxBecomesLowerBound);
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::
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

    _simplexTableau._rightHandSides[rowIdx] = SimplexTraitsT::add(
        _simplexTableau._rightHandSides[rowIdx], addedValue);
  }

  isColumnAtLowerBoundBitset[enteringColumnIdx] =
      !isColumnAtLowerBoundBitset[enteringColumnIdx];
  isColumnAtUpperBoundBitset[enteringColumnIdx] =
      !isColumnAtUpperBoundBitset[enteringColumnIdx];
}

template <typename T, typename SimplexTraitsT>
std::optional<PivotRowData<T>>
RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::chooseRowIdx(
    const int enteringColumnIdx, const std::vector<T> &enteringColumn) {
  std::optional<int> rowIdxMostRestrictiveLowerBound;
  std::optional<T> mostRestrictiveLowerBound;

  std::optional<int> rowIdxMostRestrictiveUpperBound;
  std::optional<T> mostRestrictiveUpperBound;

  for (int rowIdx = 0; rowIdx < _simplexTableau.getRowInfos().size();
       ++rowIdx) {
    if (!SimplexTraitsT::isEligibleForPivot(enteringColumn[rowIdx]))
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

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<
    T, SimplexTraitsT>::removeArtificialVariablesFromProgram() {
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
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
    _simplexTableau._constraintMatrix[rowIdx].resize(*firstArtificialIdx);

  _simplexTableau.initMatrixRepresentations();
}

template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBoundsSparse<
    T, SimplexTraitsT>::removeArtificialVariablesFromBasis() {
  std::vector<bool> shouldRowBeRemoved(_simplexTableau._rowInfos.size(), false);

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    SPDLOG_DEBUG("FOUND BASIC ARTIFICIAL ROW IDX {} COLUMN IDX {}, RHS {}",
                 rowIdx, basicVarIdx, _simplexTableau._rightHandSides[rowIdx]);

    std::optional<int> nonZeroEntryColumnIndex;
    const auto pivotRow = _simplexTableau.computeTableauRowPFISparse(rowIdx);

    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
      if (_simplexTableau._variableInfos[j]._isArtificial ||
          _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        continue;

      if (SimplexTraitsT::isEligibleForPivot(pivotRow._normalVec[j])) {
        nonZeroEntryColumnIndex = j;
        break;
      }
    }

    if (!nonZeroEntryColumnIndex.has_value()) {
      SPDLOG_INFO("ZERO ROW IDX {} - {}", rowIdx,
                  fmt::join(pivotRow._normalVec, ", "));
      shouldRowBeRemoved[rowIdx] = true;
    } else {
      _simplexTableau.pivotImplicitBoundsSparse(
          rowIdx, *nonZeroEntryColumnIndex,
          _simplexTableau.computeTableauColumnPFISparse(
              *nonZeroEntryColumnIndex),
          pivotRow, true);
    }
  }

  if (std::any_of(shouldRowBeRemoved.begin(), shouldRowBeRemoved.end(),
                  [](const bool val) { return val; })) {
    SPDLOG_INFO("REDUNDANT CONSTRAINTS IN LP FORMULATION");
    removeRows(shouldRowBeRemoved);
    return _simplexTableau.reinversionPFISparse();
  }

  return true;
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::removeRows(
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
  // FIXME PFI CASE
  for (int i = 0; i < _simplexTableau._basisMatrixInverse.size(); ++i)
    removeElements(_simplexTableau._basisMatrixInverse[i], shouldRowBeRemoved);

  removeElements(_simplexTableau._basisMatrixInverse, shouldRowBeRemoved);
}
template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<
    T, SimplexTraitsT>::lexicographicReoptimization(const bool minimize) {
  int curVarIdxToBeOptimized = 0;
  int varsFixedCount = 0;
  while (curVarIdxToBeOptimized < _simplexTableau._variableInfos.size() &&
         varsFixedCount < _simplexTableau._variableInfos.size()) {
    fixNonBasicVariables(varsFixedCount);
    if (!_simplexTableau._variableInfos[curVarIdxToBeOptimized]._isFixed) {
      SPDLOG_INFO("VAR IDX {} TO BE OPTIMIZED", curVarIdxToBeOptimized);
      _simplexTableau.setObjectiveSparse(
          singleVarObjective(curVarIdxToBeOptimized, minimize));
      runImpl();
    }
    ++curVarIdxToBeOptimized;
  }
  _simplexTableau.setObjectiveSparse(
      _simplexTableau._initialProgram.getObjective());
  unfixAllVariables();
}
template <typename T, typename SimplexTraitsT>
std::vector<T>
RevisedPrimalSimplexPFIBoundsSparse<T, SimplexTraitsT>::singleVarObjective(
    const int varIdx, const bool minimize) {
  std::vector<T> result(_simplexTableau._variableInfos.size());
  result[varIdx] = minimize ? 1.0 : -1.0;
  return result;
}
template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<
    T, SimplexTraitsT>::fixNonBasicVariables(int &varsFixedCount) {
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[varIdx] &&
        !SimplexTraitsT::equal(_simplexTableau._reducedCosts[varIdx], 0.0)) {
      _simplexTableau._variableInfos[varIdx]._isFixed = true;
      ++varsFixedCount;
    }
  }
}
template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBoundsSparse<T,
                                         SimplexTraitsT>::unfixAllVariables() {
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
        _simplexTableau._variableLowerBounds[varIdx] ==
            _simplexTableau._variableUpperBounds[varIdx])
      continue;

    _simplexTableau._variableInfos[varIdx]._isFixed = false;
  }
}

template class RevisedPrimalSimplexPFIBoundsSparse<double>;
template class RevisedPrimalSimplexPFIBoundsSparse<long double>;