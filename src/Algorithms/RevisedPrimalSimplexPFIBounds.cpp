#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"

#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexValidator.h"
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
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::RevisedPrimalSimplexPFIBounds(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const int32_t objValueLoggingFrequency, const int32_t reinversionFrequency,
    const ValidateSimplex validateSimplex)
    : _simplexTableau(simplexTableau),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency),
      _validateSimplex(validateSimplex) {}

template <typename T, typename SimplexTraitsT>
std::string RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::type() const {
  return "REVISED PRIMAL SIMPLEX (" +
         std::string(SimplexTraitsT::useSparseRepresentationValue ? "SPARSE"
                                                                  : "NORMAL") +
         ')';
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::runPhaseOne() {
  SPDLOG_INFO("BASIS SIZE {}, COLUMN PIVOT RULE {}",
              _simplexTableau._rowInfos.size(),
              primalSimplexColumnPivotRuleToStr(_primalSimplexColumnPivotRule));
  auto artLpOptStats = runImpl("PHASE_ONE");

  if (NumericalTraitsT::greater(_simplexTableau._objectiveValue, 0.0)) {
    SPDLOG_WARN(
        "PROGRAM WITH ARTIFICIAL VARIABLE HAS OPTIMUM {} GREATER THAN 0 - "
        "INITIAL PROGRAM IS INFEASIBLE",
        _simplexTableau._objectiveValue);
    artLpOptStats._phaseOneSucceeded = false;
    return artLpOptStats;
  }

  removeArtificialVariablesFromBasis();
  removeArtificialVariablesFromProgram();
  artLpOptStats._phaseOneSucceeded = _simplexTableau.reinversion();
  return artLpOptStats;
}
template <typename T, typename SimplexTraitsT>
LPOptStatistics<T>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::runPhaseTwo() {
  _simplexTableau.setObjective(_simplexTableau._initialProgram.getObjective());
  return runImpl("PHASE_TWO");
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T> RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::runImpl(
    const std::string &lpNameSuffix, const bool printSummary) {
  [[maybe_unused]] int iterCount = 1;
  SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  LPOptStatistics<T> lpOptStatistics{
      ._lpName = (_simplexTableau.getName() + '_' + lpNameSuffix),
      ._simplexAlgorithmType = type(),
      ._reinversionFrequency = _reinversionFrequency};
  while (true) {
    const bool iterResult = runOneIteration();
    if (iterResult)
      break;

    _simplexTableau.calculateCurrentObjectiveValue();
    _simplexTableau.calculateSolution();

    lpOptStatistics._consecutiveObjectiveValues.push_back(
        _simplexTableau.getCurrentObjectiveValue());

    ++iterCount;
    tryLogObjValue(iterCount);

    if (!tryValidateIteration(lpOptStatistics))
      break;

    if (!tryReinversion(iterCount, lpOptStatistics))
      break;

    if (!checkIterationLimit(iterCount))
      break;

    SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  }
  if (_simplexTableau.getLPOptResult() ==
      LPOptimizationResult::BOUNDED_AND_FEASIBLE)
    tryValidateOptimalSolutions(lpOptStatistics);

  if (printSummary)
    SPDLOG_INFO("{} ENDED, LP OPT RESULT {}, ITERATION COUNT {}", type(),
                lpOptimizationResultToStr(_simplexTableau._result), iterCount);

  lpOptStatistics._optResult = _simplexTableau.getLPOptResult();
  lpOptStatistics._optimalValue = _simplexTableau.getCurrentObjectiveValue();
  lpOptStatistics._iterationCount = iterCount;

  return lpOptStatistics;
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::tryLogObjValue(
    const int iterCount) {
  if (_objValueLoggingFrequency &&
      (iterCount % _objValueLoggingFrequency == 0)) {
    SPDLOG_INFO("ITERATION {}", iterCount);
    SPDLOG_INFO("{}\n", _simplexTableau.toStringObjectiveValue());
  }
}

template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::tryReinversion(
    const int iterCount, const LPOptStatistics<T> &lpOptStatistics) {
  if (_reinversionFrequency && (iterCount % _reinversionFrequency == 0)) {
    if (!_simplexTableau.reinversion()) {
      SPDLOG_WARN("STOPPING {} BECAUSE OF FAILED REINVERSION", type());
      _simplexTableau._result = LPOptimizationResult::FAILED_REINVERSION;
      return false;
    }
    if (!tryValidateIteration(lpOptStatistics))
      return false;
  }
  return true;
}
template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::tryValidateIteration(
    const LPOptStatistics<T> &lpOptStatistics) {
  if (_validateSimplex == ValidateSimplex::NO)
    return true;

  const auto validationResult =
      SimplexValidator<T, SimplexTraitsT>(_simplexTableau, lpOptStatistics)
          .validatePrimalIteration();
  if (!validationResult) {
    SPDLOG_ERROR("ITERATION VALIDATION FAILED - {}", validationResult.error());
    _simplexTableau._result = LPOptimizationResult::FAILED_VALIDATION;
    return false;
  }

  return true;
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::
    tryValidateOptimalSolutions(const LPOptStatistics<T> &lpOptStatistics) {
  if (_validateSimplex == ValidateSimplex::NO)
    return;

  const auto validationResult =
      SimplexValidator<T, SimplexTraitsT>(_simplexTableau, lpOptStatistics)
          .validateOptimality(SimplexType::PRIMAL);
  if (!validationResult) {
    SPDLOG_ERROR("OPTIMALITY VALIDATION FAILED - {}", validationResult.error());
    _simplexTableau._result = LPOptimizationResult::FAILED_VALIDATION;
  }
}

template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::checkIterationLimit(
    const int iterCount) {
  constexpr size_t HARD_ITERATION_LIMIT = 100000;
  if (iterCount > HARD_ITERATION_LIMIT) {
    _simplexTableau._result = LPOptimizationResult::REACHED_ITERATION_LIMIT;
    return false;
  }
  return true;
}

template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::runOneIteration() {
  const std::optional<int> enteringColumnIdx = chooseEnteringColumn();
  if (!enteringColumnIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

  SPDLOG_DEBUG("ENTERING COLUMN IDX {} REDUCED COST {}", *enteringColumnIdx,
               _simplexTableau._reducedCosts[*enteringColumnIdx]);

  const auto enteringColumn =
      _simplexTableau.computeTableauColumnGeneric(*enteringColumnIdx);
  const std::optional<PivotRowData<T>> pivotRowData =
      chooseRowIdx(*enteringColumnIdx, enteringColumn);
  if (!pivotRowData.has_value()) {
    _simplexTableau._result = LPOptimizationResult::UNBOUNDED;
    return true;
  }

  changeTableau(*pivotRowData, *enteringColumnIdx, enteringColumn);
  return false;
}

template <typename T, typename SimplexTraitsT>
std::optional<int>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::chooseEnteringColumn() {
  switch (_primalSimplexColumnPivotRule) {
  case PrimalSimplexColumnPivotRule::FIRST_ELIGIBLE:
    return chooseEnteringColumnFirstEligible();
  case PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST:
    return chooseEnteringColumnBiggestAbsReducedCost();
  }
}

template <typename T, typename SimplexTraitsT>
std::optional<int> RevisedPrimalSimplexPFIBounds<
    T, SimplexTraitsT>::chooseEnteringColumnFirstEligible() {
  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size();
       ++columnIdx) {
    SPDLOG_TRACE(
        "COL IDX {} ART {} BASIC {} RED COST {} < 0.0 {} > 0.0 {}", columnIdx,
        _simplexTableau._variableInfos[columnIdx]._isArtificial,
        _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx],
        _simplexTableau._reducedCosts[columnIdx],
        NumericalTraitsT::less(_simplexTableau._reducedCosts[columnIdx], 0.0),
        NumericalTraitsT::greater(_simplexTableau._reducedCosts[columnIdx],
                                  0.0));

    SPDLOG_TRACE("COL IDX {}, IS LB {}, IS UB {}", columnIdx,
                 _simplexTableau._simplexBasisData
                     ._isColumnAtLowerBoundBitset[columnIdx],
                 _simplexTableau._simplexBasisData
                     ._isColumnAtUpperBoundBitset[columnIdx]);

    if (!_simplexTableau.isColumnAllowedToEnterBasis(columnIdx))
      continue;

    if (_simplexTableau._reducedCosts[columnIdx] <
            -NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE &&
        _simplexTableau._simplexBasisData
            ._isColumnAtLowerBoundBitset[columnIdx])
      return columnIdx;

    if (_simplexTableau._reducedCosts[columnIdx] >
            NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE &&
        _simplexTableau._simplexBasisData
            ._isColumnAtUpperBoundBitset[columnIdx])
      return columnIdx;
  }
  return std::nullopt;
}

template <typename T, typename SimplexTraitsT>
std::optional<int> RevisedPrimalSimplexPFIBounds<
    T, SimplexTraitsT>::chooseEnteringColumnBiggestAbsReducedCost() {
  std::optional<int> bestColumnIdx;
  std::optional<T> biggestAbsReducedCost;

  const auto tryUpdateBest = [&](const int colIdx) {
    if ((_simplexTableau._simplexBasisData
             ._isColumnAtLowerBoundBitset[colIdx] &&
         _simplexTableau._reducedCosts[colIdx] <
             -NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE) ||
        (_simplexTableau._simplexBasisData
             ._isColumnAtUpperBoundBitset[colIdx] &&
         _simplexTableau._reducedCosts[colIdx] >
             NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE)) {
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
void RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::changeTableau(
    const PivotRowData<T> &pivotRowData, const int enteringColumnIdx,
    const VectorT &enteringColumn) {
  if (pivotRowData._noBasisChangeNeeded)
    moveVarFromOneBoundToAnother(pivotRowData, enteringColumnIdx,
                                 enteringColumn);
  else
    _simplexTableau.pivotImplicitBoundsGeneric(
        *pivotRowData._pivotRowIdx, enteringColumnIdx, enteringColumn,
        _simplexTableau.computeTableauRowGeneric(*pivotRowData._pivotRowIdx),
        pivotRowData._departingIdxBecomesLowerBound);
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::
    moveVarFromOneBoundToAnother(const PivotRowData<T> &pivotRowData,
                                 const int enteringColumnIdx,
                                 const VectorT &enteringColumn) {
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

    _simplexTableau._rightHandSides[rowIdx] = NumericalTraitsT::add(
        _simplexTableau._rightHandSides[rowIdx], addedValue);
  }

  isColumnAtLowerBoundBitset[enteringColumnIdx] =
      !isColumnAtLowerBoundBitset[enteringColumnIdx];
  isColumnAtUpperBoundBitset[enteringColumnIdx] =
      !isColumnAtUpperBoundBitset[enteringColumnIdx];
}

template <typename T, typename SimplexTraitsT>
std::optional<PivotRowData<T>>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::chooseRowIdx(
    const int enteringColumnIdx, const VectorT &enteringColumn) {
  std::optional<int> rowIdxMostRestrictiveLowerBound;
  std::optional<T> mostRestrictiveLowerBound;

  std::optional<int> rowIdxMostRestrictiveUpperBound;
  std::optional<T> mostRestrictiveUpperBound;

  for (int rowIdx = 0; rowIdx < _simplexTableau.getRowInfos().size();
       ++rowIdx) {
    if (!NumericalTraitsT::isEligibleForPivot(enteringColumn[rowIdx]))
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
void RevisedPrimalSimplexPFIBounds<
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
  _simplexTableau._variableLowerBounds.resize(*firstArtificialIdx);
  _simplexTableau._variableUpperBounds.resize(*firstArtificialIdx);
  _simplexTableau._simplexBasisData.resizeVarCount(*firstArtificialIdx);
  _simplexTableau._reducedCosts.resize(*firstArtificialIdx);
  _simplexTableau._x.resize(*firstArtificialIdx);

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
    _simplexTableau._constraintMatrix[rowIdx].resize(*firstArtificialIdx);

  _simplexTableau.initMatrixRepresentations();
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBounds<
    T, SimplexTraitsT>::removeArtificialVariablesFromBasis() {
  std::vector<bool> shouldRowBeRemoved(_simplexTableau._rowInfos.size(), false);

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    SPDLOG_DEBUG("FOUND BASIC ARTIFICIAL ROW IDX {} COLUMN IDX {}, RHS {}",
                 rowIdx, basicVarIdx, _simplexTableau._rightHandSides[rowIdx]);

    std::optional<int> nonZeroEntryColumnIndex;
    const auto pivotRow = _simplexTableau.computeTableauRowGeneric(rowIdx);

    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
      if (_simplexTableau._variableInfos[j]._isArtificial ||
          _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        continue;

      if (NumericalTraitsT::isEligibleForPivot(pivotRow[j])) {
        nonZeroEntryColumnIndex = j;
        break;
      }
    }

    if (!nonZeroEntryColumnIndex.has_value()) {
      shouldRowBeRemoved[rowIdx] = true;
    } else {
      _simplexTableau.pivotImplicitBoundsGeneric(
          rowIdx, *nonZeroEntryColumnIndex,
          _simplexTableau.computeTableauColumnGeneric(*nonZeroEntryColumnIndex),
          pivotRow, true);
    }
  }

  removeRows(shouldRowBeRemoved);
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::removeRows(
    const std::vector<bool> &shouldRowBeRemoved) {
  const auto rowsToBeRemoved =
      std::count(shouldRowBeRemoved.begin(), shouldRowBeRemoved.end(), true);
  if (rowsToBeRemoved == 0)
    return;

  SPDLOG_INFO("REDUNDANT {} CONSTRAINTS IN LP FORMULATION", rowsToBeRemoved);

  auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset,
         isColumnAtLowerBoundBitset, _2] = _simplexTableau._simplexBasisData;

  std::vector<int> oldRowIdxToNewRowIdx(shouldRowBeRemoved.size());
  int curNewRowIdx = 0;

  for (int rowIdx = 0; rowIdx < shouldRowBeRemoved.size(); ++rowIdx)
    if (shouldRowBeRemoved[rowIdx]) {
      const auto basicColumnIdx = rowToBasisColumnIdxMap[rowIdx];
      isBasicColumnIndexBitset[basicColumnIdx] = false;
      isColumnAtLowerBoundBitset[basicColumnIdx] = true;
    } else
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

template class RevisedPrimalSimplexPFIBounds<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class RevisedPrimalSimplexPFIBounds<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class RevisedPrimalSimplexPFIBounds<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class RevisedPrimalSimplexPFIBounds<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
