#include "src/Algorithms/DualSimplex.h"

#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexValidator.h"
#include "src/Util/SpdlogHeader.h"
#include "src/Util/Time.h"

#include <fmt/format.h>

template <typename T, typename SimplexTraitsT>
DualSimplex<T, SimplexTraitsT>::DualSimplex(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const DualSimplexRowPivotRule dualSimplexRowPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager),
      _dualSimplexRowPivotRule(dualSimplexRowPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption) {}

template <typename T, typename SimplexTraitsT>
std::string DualSimplex<T, SimplexTraitsT>::type() const {
  return fmt::format(
      "DUAL SIMPLEX ({}, {})",
      simplexTableauTypeToStr(_simplexTableau._simplexTableauType),
      matrixRepresentationTypeToStr(SimplexTraitsT::matrixRepresentationType));
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T> DualSimplex<T, SimplexTraitsT>::run(
    const std::string &lpNameSuffix,
    const PrintSimplexOptSummary printSimplexOptSummary,
    const DualPhase dualPhase) {
  if (printSimplexOptSummary == PrintSimplexOptSummary::YES) {
    SPDLOG_INFO("LP NAME {} BASIS SIZE {}, ROW PIVOT RULE {}",
                _simplexTableau._initialProgram.getName(),
                _simplexTableau._rowInfos.size(),
                dualSimplexRowPivotRuleToStr(_dualSimplexRowPivotRule));
  }
  SPDLOG_TRACE("{}\n", _simplexTableau.toString());

  LPOptStatistics<T> lpOptStatistics{
      ._lpName = _simplexTableau.getName() +
                 (!lpNameSuffix.empty() ? "_" + lpNameSuffix : ""),
      ._algorithmType = type(),
      ._reinversionFrequency = _reinversionManager.reinversionFrequency()};
  [[maybe_unused]] int iterCount = 1;
  lpOptStatistics._elapsedTimeSec = executeAndMeasureTime([&] {
    while (true) {
      const bool iterResult = runOneIteration();
      if (iterResult)
        break;

      _simplexTableau.calculateSolution();
      _simplexTableau.calculateCurrentObjectiveValue();

      lpOptStatistics._consecutiveObjectiveValues.push_back(
          _simplexTableau._objectiveValue);

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

  if (_simplexTableau._result == LPOptimizationResult::BOUNDED_AND_FEASIBLE)
    tryValidateOptimalSolutions(lpOptStatistics, dualPhase);

  if (printSimplexOptSummary == PrintSimplexOptSummary::YES) {
    SPDLOG_INFO("{} ENDED", type());
    SPDLOG_INFO("LP OPT RESULT {}, OPT VALUE {}",
                lpOptimizationResultToStr(_simplexTableau._result),
                _simplexTableau._objectiveValue);
    SPDLOG_INFO("ELAPSED TIME {} SECONDS, ITERATION COUNT {}",
                lpOptStatistics._elapsedTimeSec, iterCount);
  }

  lpOptStatistics._optResult = _simplexTableau._result;
  if (lpOptStatistics._optResult == LPOptimizationResult::FAILED_VALIDATION) {
    SPDLOG_ERROR("OPTIMALITY VALIDATION FAILED");
  }
  lpOptStatistics._optimalValue = _simplexTableau._objectiveValue;
  lpOptStatistics._iterationCount = iterCount;

  return lpOptStatistics;
}

template <typename T, typename SimplexTraitsT>
void DualSimplex<T, SimplexTraitsT>::tryLogObjValue(const int iterCount) {
  if (_objValueLoggingFrequency &&
      (iterCount % _objValueLoggingFrequency == 0)) {
    SPDLOG_INFO("ITERATION {}", iterCount);
    SPDLOG_INFO("{}\n", _simplexTableau.toStringObjectiveValue());
  }
}

template <typename T, typename SimplexTraitsT>
bool DualSimplex<T, SimplexTraitsT>::tryReinversion(
    const int iterCount, const LPOptStatistics<T> &lpOptStatistics) {
  if (!_reinversionManager.tryReinverse()) {
    SPDLOG_WARN("STOPPING {} BECAUSE OF FAILED REINVERSION", type());
    _simplexTableau._result = LPOptimizationResult::FAILED_REINVERSION;
    return false;
  }

  return tryValidateIteration(iterCount, lpOptStatistics);
}

template <typename T, typename SimplexTraitsT>
bool DualSimplex<T, SimplexTraitsT>::tryValidateIteration(
    const int iterCount, const LPOptStatistics<T> &lpOptStatistics) {
  if (_validateSimplexOption == ValidateSimplexOption::DONT_VALIDATE)
    return true;

  const auto validationResult =
      SimplexValidator<T, SimplexTraitsT>(_simplexTableau, lpOptStatistics)
          .validateDualIteration();
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
void DualSimplex<T, SimplexTraitsT>::tryValidateOptimalSolutions(
    const LPOptStatistics<T> &lpOptStatistics, const DualPhase dualPhase) {
  if (_validateSimplexOption == ValidateSimplexOption::DONT_VALIDATE)
    return;

  const auto validationResult =
      SimplexValidator<T, SimplexTraitsT>(_simplexTableau, lpOptStatistics)
          .validateOptimality(SimplexType::DUAL, dualPhase);
  if (!validationResult) {
    SPDLOG_ERROR("OPTIMALITY VALIDATION FAILED - {}", validationResult.error());
    if (_validateSimplexOption ==
        ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR) {
      _simplexTableau._result = LPOptimizationResult::FAILED_VALIDATION;
    }
  }
}

template <typename T, typename SimplexTraitsT>
bool DualSimplex<T, SimplexTraitsT>::checkIterationLimit(const int iterCount) {
  constexpr size_t HARD_ITERATION_LIMIT = 500000;
  if (iterCount > HARD_ITERATION_LIMIT) {
    _simplexTableau._result = LPOptimizationResult::REACHED_ITERATION_LIMIT;
    return false;
  }
  return true;
}

template <typename T, typename SimplexTraitsT>
bool DualSimplex<T, SimplexTraitsT>::checkObjectiveProgress(
    const LPOptStatistics<T> &lpOptStatistics) {
  constexpr size_t ITERATION_WINDOW_SIZE = 30000;
  if (lpOptStatistics._consecutiveObjectiveValues.size() >=
      ITERATION_WINDOW_SIZE) {
    const auto currentObjValue =
        lpOptStatistics._consecutiveObjectiveValues.back();
    const auto objValueBeforeTheWindow =
        lpOptStatistics._consecutiveObjectiveValues
            [lpOptStatistics._consecutiveObjectiveValues.size() -
             ITERATION_WINDOW_SIZE];
    if (currentObjValue <=
        objValueBeforeTheWindow +
            NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE) {
      _simplexTableau._result = LPOptimizationResult::TOO_SMALL_PROGRESS;
      return false;
    }
  }
  return true;
}

template <typename T, typename SimplexTraitsT>
bool DualSimplex<T, SimplexTraitsT>::runOneIteration() {
  const std::optional<int> pivotRowIdx = chooseRow();
  if (!pivotRowIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

  SPDLOG_DEBUG("PIVOT ROW IDX {} RHS VALUE {}", *pivotRowIdx,
               _simplexTableau._rightHandSides[*pivotRowIdx]);
  const auto basicColumnIdx =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[*pivotRowIdx];
  const bool isPivotRowUnderLowerBound =
      _simplexTableau._rightHandSides[*pivotRowIdx] <
      *_simplexTableau._variableLowerBounds[basicColumnIdx] -
          NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE;
  return pivot(*pivotRowIdx, std::nullopt, isPivotRowUnderLowerBound);
}

template <typename T, typename SimplexTraitsT>
bool DualSimplex<T, SimplexTraitsT>::pivot(
    const int pivotRowIdx, const std::optional<int> customEnteringColumnIdx,
    const bool isPivotRowUnderLowerBound) {
  const auto pivotRowSharedPtr =
      _simplexTableau.computeTableauRowGeneric(pivotRowIdx);
  const auto &pivotRow = *pivotRowSharedPtr;
  const std::optional<int> enteringColumnIdx =
      customEnteringColumnIdx.has_value()
          ? customEnteringColumnIdx
          : chooseEnteringColumnIdx(pivotRowIdx, pivotRow,
                                    isPivotRowUnderLowerBound);
  if (!enteringColumnIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::INFEASIBLE;
    return true;
  }

  const auto enteringColumn =
      _simplexTableau.computeTableauColumnGeneric(*enteringColumnIdx);

  _simplexTableau.pivotImplicitBoundsGeneric(pivotRowIdx, *enteringColumnIdx,
                                             enteringColumn, pivotRow,
                                             isPivotRowUnderLowerBound);
  return false;
}

template <typename T, typename SimplexTraitsT>
std::optional<int> DualSimplex<T, SimplexTraitsT>::chooseRow() {
  switch (_dualSimplexRowPivotRule) {
  case DualSimplexRowPivotRule::FIRST_ELIGIBLE:
    return chooseRowFirstEligible();
  case DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION:
    return chooseRowBiggestViolation();
  }
}

template <typename T, typename SimplexTraitsT>
std::optional<int> DualSimplex<T, SimplexTraitsT>::chooseRowFirstEligible() {
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicColumnIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    const auto lowerBound =
        _simplexTableau._variableLowerBounds[basicColumnIdx];
    if (lowerBound.has_value() &&
        _simplexTableau._rightHandSides[rowIdx] <
            *lowerBound - NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)
      return rowIdx;

    const auto upperBound =
        _simplexTableau._variableUpperBounds[basicColumnIdx];
    if (upperBound.has_value() &&
        _simplexTableau._rightHandSides[rowIdx] >
            *upperBound + NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)
      return rowIdx;
  }
  return std::nullopt;
}

template <typename T, typename SimplexTraitsT>
std::optional<int> DualSimplex<T, SimplexTraitsT>::chooseRowBiggestViolation() {
  std::optional<int> bestRowIdx;
  std::optional<T> biggestViolation;

  const auto tryUpdateBest = [&](const int rowIdx, const T curViolation) {
    if (!biggestViolation.has_value() || (*biggestViolation < curViolation)) {
      bestRowIdx = rowIdx;
      biggestViolation = curViolation;
    }
  };

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicColumnIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    const auto lowerBound =
        _simplexTableau._variableLowerBounds[basicColumnIdx];
    if (lowerBound.has_value() &&
        _simplexTableau._rightHandSides[rowIdx] <
            *lowerBound - NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)
      tryUpdateBest(rowIdx,
                    *lowerBound - _simplexTableau._rightHandSides[rowIdx]);

    const auto upperBound =
        _simplexTableau._variableUpperBounds[basicColumnIdx];
    if (upperBound.has_value() &&
        _simplexTableau._rightHandSides[rowIdx] >
            *upperBound + NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)
      tryUpdateBest(rowIdx,
                    _simplexTableau._rightHandSides[rowIdx] - *upperBound);
  }

  return bestRowIdx;
}

template <typename T, typename SimplexTraitsT>
std::optional<int> DualSimplex<T, SimplexTraitsT>::chooseEnteringColumnIdx(
    const int pivotRowIdx, const VectorT &pivotRow,
    const bool isPivotRowUnderLowerBound) {
  std::optional<int> mostRestrictiveColumnIdx;
  std::optional<T> mostRestrictiveColumnBound;

  auto &isColumnAtLowerBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset;

  const auto currentRestrictionBound =
      [&](const int colIdx) -> std::optional<T> {
    if (NumericalTraitsT::isEligibleForPivot(pivotRow[colIdx])) {
      if (isPivotRowUnderLowerBound) {
        if (isColumnAtLowerBoundBitset[colIdx]) {
          if (NumericalTraitsT::less(pivotRow[colIdx], 0.0))
            return _simplexTableau._reducedCosts[colIdx] / (-pivotRow[colIdx]);
        } else {
          if (NumericalTraitsT::greater(pivotRow[colIdx], 0.0))
            return (-_simplexTableau._reducedCosts[colIdx]) / pivotRow[colIdx];
        }
      } else {
        if (isColumnAtLowerBoundBitset[colIdx]) {
          if (NumericalTraitsT::greater(pivotRow[colIdx], 0.0))
            return _simplexTableau._reducedCosts[colIdx] / pivotRow[colIdx];
        } else {
          if (NumericalTraitsT::less(pivotRow[colIdx], 0.0))
            return (-_simplexTableau._reducedCosts[colIdx]) /
                   (-pivotRow[colIdx]);
        }
      }
    }

    return std::nullopt;
  };

  const auto tryUpdateBest = [&](const int colIdx) {
    if (const auto curBound = currentRestrictionBound(colIdx);
        curBound.has_value()) {
      if (!mostRestrictiveColumnBound.has_value() ||
          (*curBound < *mostRestrictiveColumnBound)) {
        mostRestrictiveColumnIdx = colIdx;
        mostRestrictiveColumnBound = *curBound;
      }
    }
  };

  for (int columnIdx = 0; columnIdx < _simplexTableau._variableInfos.size();
       ++columnIdx)
    if (_simplexTableau.isColumnAllowedToEnterBasis(columnIdx))
      tryUpdateBest(columnIdx);

  return mostRestrictiveColumnIdx;
}

template class DualSimplex<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplex<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class DualSimplex<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class DualSimplex<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
