#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"

#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexValidator.h"
#include "src/Util/SpdlogHeader.h"

template <typename T, typename SimplexTraitsT>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::RevisedDualSimplexPFIBounds(
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
std::string RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::type() const {
  return "REVISED DUAL SIMPLEX (" +
         std::string(SimplexTraitsT::useSparseRepresentationValue ? "SPARSE"
                                                                  : "NORMAL") +
         ')';
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T> RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::run(
    const std::string &lpNameSuffix) {
  SPDLOG_INFO("BASIS SIZE {} COLUMN PIVOT RULE {}",
              _simplexTableau._rowInfos.size(),
              dualSimplexRowPivotRuleToStr(_dualSimplexRowPivotRule));
  SPDLOG_TRACE("{}\n", _simplexTableau.toString());

  LPOptStatistics<T> lpOptStatistics{
      ._lpName = _simplexTableau.getName() + '_' + lpNameSuffix,
      ._simplexAlgorithmType = type(),
      ._reinversionFrequency = _reinversionManager.reinversionFrequency()};
  [[maybe_unused]] int iterCount = 1;
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

    if (!checkIterationLimit(iterCount))
      break;

    SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  }

  if (_simplexTableau.getLPOptResult() ==
      LPOptimizationResult::BOUNDED_AND_FEASIBLE)
    tryValidateOptimalSolutions(lpOptStatistics);

  SPDLOG_INFO("{} ENDED, LP OPT RESULT {}, ITERATION COUNT {}", type(),
              lpOptimizationResultToStr(_simplexTableau._result), iterCount);

  lpOptStatistics._optResult = _simplexTableau.getLPOptResult();
  lpOptStatistics._optimalValue = _simplexTableau.getCurrentObjectiveValue();
  lpOptStatistics._iterationCount = iterCount;

  return lpOptStatistics;
}

template <typename T, typename SimplexTraitsT>
void RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::tryLogObjValue(
    const int iterCount) {
  if (_objValueLoggingFrequency &&
      (iterCount % _objValueLoggingFrequency == 0)) {
    SPDLOG_INFO("ITERATION {}", iterCount);
    SPDLOG_INFO("{}\n", _simplexTableau.toStringObjectiveValue());
  }
}

template <typename T, typename SimplexTraitsT>
bool RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::tryReinversion(
    const int iterCount, const LPOptStatistics<T> &lpOptStatistics) {
  if (!_reinversionManager.tryReinverse()) {
    SPDLOG_WARN("STOPPING {} BECAUSE OF FAILED REINVERSION", type());
    _simplexTableau._result = LPOptimizationResult::FAILED_REINVERSION;
    return false;
  }

  return tryValidateIteration(iterCount, lpOptStatistics);
}

template <typename T, typename SimplexTraitsT>
bool RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::tryValidateIteration(
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
void RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::
    tryValidateOptimalSolutions(const LPOptStatistics<T> &lpOptStatistics) {
  if (_validateSimplexOption == ValidateSimplexOption::DONT_VALIDATE)
    return;

  const auto validationResult =
      SimplexValidator<T, SimplexTraitsT>(_simplexTableau, lpOptStatistics)
          .validateOptimality(SimplexType::DUAL);
  if (!validationResult) {
    SPDLOG_ERROR("OPTIMALITY VALIDATION FAILED - {}", validationResult.error());
    if (_validateSimplexOption ==
        ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR) {
      _simplexTableau._result = LPOptimizationResult::FAILED_VALIDATION;
    }
  }
}

template <typename T, typename SimplexTraitsT>
bool RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::checkIterationLimit(
    const int iterCount) {
  constexpr size_t HARD_ITERATION_LIMIT = 500000;
  if (iterCount > HARD_ITERATION_LIMIT) {
    _simplexTableau._result = LPOptimizationResult::REACHED_ITERATION_LIMIT;
    return false;
  }
  return true;
}

template <typename T, typename SimplexTraitsT>
bool RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::runOneIteration() {
  const std::optional<int> pivotRowIdx = chooseRow();
  if (!pivotRowIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

  SPDLOG_DEBUG("PIVOT ROW IDX {} RHS VALUE {}", *pivotRowIdx,
               _simplexTableau._rightHandSides[*pivotRowIdx]);
  const auto pivotRow = _simplexTableau.computeTableauRowGeneric(*pivotRowIdx);
  const auto basicColumnIdx =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[*pivotRowIdx];
  const bool isPivotRowUnderLowerBound =
      _simplexTableau._rightHandSides[*pivotRowIdx] <
      *_simplexTableau._variableLowerBounds[basicColumnIdx] -
          NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE;
  const auto enteringColumnIdx = chooseEnteringColumnIdx(
      *pivotRowIdx, pivotRow, isPivotRowUnderLowerBound);
  if (!enteringColumnIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::INFEASIBLE;
    return true;
  }

  const auto enteringColumn =
      _simplexTableau.computeTableauColumnGeneric(*enteringColumnIdx);

  _simplexTableau.pivotImplicitBoundsGeneric(*pivotRowIdx, *enteringColumnIdx,
                                             enteringColumn, pivotRow,
                                             isPivotRowUnderLowerBound);
  return false;
}

template <typename T, typename SimplexTraitsT>
std::optional<int> RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::chooseRow() {
  switch (_dualSimplexRowPivotRule) {
  case DualSimplexRowPivotRule::FIRST_ELIGIBLE:
    return chooseRowFirstEligible();
  case DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION:
    return chooseRowBiggestViolation();
  }
}

template <typename T, typename SimplexTraitsT>
std::optional<int>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::chooseRowFirstEligible() {
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
std::optional<int>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::chooseRowBiggestViolation() {
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
std::optional<int>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::chooseEnteringColumnIdx(
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

  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size();
       ++columnIdx)
    if (_simplexTableau.isColumnAllowedToEnterBasis(columnIdx))
      tryUpdateBest(columnIdx);

  return mostRestrictiveColumnIdx;
}

template class RevisedDualSimplexPFIBounds<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class RevisedDualSimplexPFIBounds<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class RevisedDualSimplexPFIBounds<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class RevisedDualSimplexPFIBounds<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
