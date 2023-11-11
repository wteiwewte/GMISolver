#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"

#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexTableauResizer.h"
#include "src/Algorithms/SimplexValidator.h"
#include "src/DataModel/LinearProgram.h"
#include "src/Util/SpdlogHeader.h"
#include "src/Util/Time.h"

template <typename T, typename SimplexTraitsT>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::RevisedPrimalSimplexPFIBounds(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption) {}

template <typename T, typename SimplexTraitsT>
std::string RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::type() const {
  return fmt::format(
      "REVISED PRIMAL SIMPLEX ({}, {})",
      simplexTableauTypeToStr(_simplexTableau._simplexTableauType),
      matrixRepresentationTypeToStr(SimplexTraitsT::matrixRepresentationType));
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::runPhaseOne() {
  SPDLOG_INFO("LP NAME {} BASIS SIZE {}, COLUMN PIVOT RULE {}",
              _simplexTableau._initialProgram.getName(),
              _simplexTableau._rowInfos.size(),
              primalSimplexColumnPivotRuleToStr(_primalSimplexColumnPivotRule));
  auto artLpOptStats = runImpl("PHASE_ONE");

  if (_simplexTableau._objectiveValue >
      NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE) {
    SPDLOG_WARN(
        "PROGRAM WITH ARTIFICIAL VARIABLE HAS OPTIMUM {} GREATER THAN 0 - "
        "INITIAL PROGRAM IS INFEASIBLE",
        _simplexTableau._objectiveValue);
    artLpOptStats._phaseOneSucceeded = false;
    return artLpOptStats;
  }

  SimplexTableauResizer simplexTableauResizer(_simplexTableau,
                                              _reinversionManager);
  artLpOptStats._phaseOneSucceeded =
      simplexTableauResizer.removeArtificialVariablesFromProgram();
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
      ._reinversionFrequency = _reinversionManager.reinversionFrequency()};
  lpOptStatistics._elapsedTimeSec = executeAndMeasureTime([&] {
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

      if (!tryValidateIteration(iterCount, lpOptStatistics))
        break;

      if (!tryReinversion(iterCount, lpOptStatistics))
        break;

      if (!checkObjectiveProgress(lpOptStatistics))
        break;

      if (!checkIterationLimit(iterCount))
        break;

      SPDLOG_TRACE("{}\n", _simplexTableau.toString());
    }
  });

  if (_simplexTableau.getLPOptResult() ==
      LPOptimizationResult::BOUNDED_AND_FEASIBLE)
    tryValidateOptimalSolutions(lpOptStatistics);

  if (printSummary) {
    SPDLOG_INFO("{} ENDED", type());
    SPDLOG_INFO("LP OPT RESULT {}, OPT VALUE {}",
                lpOptimizationResultToStr(_simplexTableau._result),
                _simplexTableau.getCurrentObjectiveValue());
    SPDLOG_INFO("ELAPSED TIME {} SECONDS, ITERATION COUNT {}",
                lpOptStatistics._elapsedTimeSec, iterCount);
  }

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
  if (!_reinversionManager.tryReinverse()) {
    SPDLOG_WARN("STOPPING {} BECAUSE OF FAILED REINVERSION", type());
    _simplexTableau._result = LPOptimizationResult::FAILED_REINVERSION;
    return false;
  }

  return tryValidateIteration(iterCount, lpOptStatistics);
}
template <typename T, typename SimplexTraitsT>
bool RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::tryValidateIteration(
    const int iterCount, const LPOptStatistics<T> &lpOptStatistics) {
  if (_validateSimplexOption == ValidateSimplexOption::DONT_VALIDATE)
    return true;

  const auto validationResult =
      SimplexValidator<T, SimplexTraitsT>(_simplexTableau, lpOptStatistics)
          .validatePrimalIteration();
  if (!validationResult) {
    SPDLOG_ERROR("ITERATION {} VALIDATION FAILED - {}", iterCount,
                 validationResult.error());

    if (_validateSimplexOption ==
        ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR) {
      _simplexTableau._result = LPOptimizationResult::FAILED_VALIDATION;
      return false;
    }
  }

  return true;
}

template <typename T, typename SimplexTraitsT>
void RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::
    tryValidateOptimalSolutions(const LPOptStatistics<T> &lpOptStatistics) {
  if (_validateSimplexOption == ValidateSimplexOption::DONT_VALIDATE)
    return;

  const auto validationResult =
      SimplexValidator<T, SimplexTraitsT>(_simplexTableau, lpOptStatistics)
          .validateOptimality(SimplexType::PRIMAL);
  if (!validationResult) {
    SPDLOG_ERROR("OPTIMALITY VALIDATION FAILED - {}", validationResult.error());
    if (_validateSimplexOption ==
        ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR) {
      _simplexTableau._result = LPOptimizationResult::FAILED_VALIDATION;
    }
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
bool RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>::checkObjectiveProgress(
    const LPOptStatistics<T> &lpOptStatistics) {
  constexpr size_t ITERATION_WINDOW_SIZE = 5000;
  if (lpOptStatistics._consecutiveObjectiveValues.size() >=
      ITERATION_WINDOW_SIZE) {
    const auto currentObjValue =
        lpOptStatistics._consecutiveObjectiveValues.back();
    const auto objValueBeforeTheWindow =
        lpOptStatistics._consecutiveObjectiveValues
            [lpOptStatistics._consecutiveObjectiveValues.size() -
             ITERATION_WINDOW_SIZE];
    if (currentObjValue >=
        objValueBeforeTheWindow -
            NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE) {
      _simplexTableau._result = LPOptimizationResult::TOO_SMALL_PROGRESS;
      return false;
    }
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

    if (std::fabs(_simplexTableau._reducedCosts[columnIdx]) >
            NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE &&
        _simplexTableau._variableInfos[columnIdx]._isFree)
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
             NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE) ||
        (std::fabs(_simplexTableau._reducedCosts[colIdx]) >
             NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE &&
         _simplexTableau._variableInfos[colIdx]._isFree)) {
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
  else {
    const auto pivotRowSharedPtr =
        _simplexTableau.computeTableauRowGeneric(*pivotRowData._pivotRowIdx);
    const auto &pivotRow = *pivotRowSharedPtr;
    _simplexTableau.pivotImplicitBoundsGeneric(
        *pivotRowData._pivotRowIdx, enteringColumnIdx, enteringColumn, pivotRow,
        pivotRowData._departingIdxBecomesLowerBound);
  }
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
  if (_simplexTableau._isVariableFreeBitset[enteringColumnIdx]) {
    SPDLOG_ERROR("FREE VAR {} SHOULDN'T BE MOVED FROM ONE BOUND TO ANOTHER",
                 enteringColumnIdx);
  }
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

    if (!lowerBound.has_value() && !upperBound.has_value())
      continue;

    if (const bool isEnteringVarIncreasing =
            _simplexTableau._reducedCosts[enteringColumnIdx] < 0.0;
        isEnteringVarIncreasing) {
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

template class RevisedPrimalSimplexPFIBounds<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class RevisedPrimalSimplexPFIBounds<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class RevisedPrimalSimplexPFIBounds<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class RevisedPrimalSimplexPFIBounds<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
