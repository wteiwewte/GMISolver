#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"

#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/SpdlogHeader.h"

template <typename T, typename SimplexTraitsT>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::RevisedDualSimplexPFIBounds(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    const DualSimplexRowPivotRule dualSimplexRowPivotRule,
    const int32_t objValueLoggingFrequency, const int32_t reinversionFrequency)
    : _simplexTableau(simplexTableau),
      _dualSimplexRowPivotRule(dualSimplexRowPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency) {
  _simplexTableau.calculateRHS();
  _simplexTableau.calculateCurrentObjectiveValue();
  _simplexTableau.calculateSolution();
}

template <typename T, typename SimplexTraitsT>
void RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::run() {
  SPDLOG_INFO("BASIS SIZE {} COLUMN PIVOT RULE {}",
              _simplexTableau._rowInfos.size(),
              dualSimplexRowPivotRuleToStr(_dualSimplexRowPivotRule));
  SPDLOG_TRACE("{}\n", _simplexTableau.toString());

  [[maybe_unused]] int iterCount = 1;
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

    SPDLOG_TRACE("{}\n", _simplexTableau.toString());
  }
  SPDLOG_INFO("DUAL SIMPLEX ENDED, ITERATION COUNT {}", iterCount);
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
  const std::vector<T> pivotRow =
      _simplexTableau.computeTableauRow(*pivotRowIdx);
  const auto basicColumnIdx =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[*pivotRowIdx];
  const bool isPivotRowUnderLowerBound = SimplexTraitsT::less(
      _simplexTableau._rightHandSides[*pivotRowIdx],
      *_simplexTableau._variableLowerBounds[basicColumnIdx]);
  const auto enteringColumnIdx = chooseEnteringColumnIdx(
      *pivotRowIdx, pivotRow, isPivotRowUnderLowerBound);
  if (!enteringColumnIdx.has_value()) {
    _simplexTableau._result = LPOptimizationResult::INFEASIBLE;
    return true;
  }

  const std::vector<T> enteringColumn =
      _simplexTableau.computeTableauColumn(*enteringColumnIdx);

  _simplexTableau.pivotImplicitBounds(*pivotRowIdx, *enteringColumnIdx,
                                      enteringColumn, pivotRow,
                                      isPivotRowUnderLowerBound);
  return false;
}

template <typename T, typename SimplexTraitsT>
std::optional<int>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::chooseRow() {
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
        SimplexTraitsT::less(_simplexTableau._rightHandSides[rowIdx],
                                *lowerBound))
      return rowIdx;

    const auto upperBound =
        _simplexTableau._variableUpperBounds[basicColumnIdx];
    if (upperBound.has_value() &&
        SimplexTraitsT::greater(_simplexTableau._rightHandSides[rowIdx],
                                   *upperBound))
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
        SimplexTraitsT::less(_simplexTableau._rightHandSides[rowIdx],
                                *lowerBound))
      tryUpdateBest(rowIdx,
                    *lowerBound - _simplexTableau._rightHandSides[rowIdx]);

    const auto upperBound =
        _simplexTableau._variableUpperBounds[basicColumnIdx];
    if (upperBound.has_value() &&
        SimplexTraitsT::greater(_simplexTableau._rightHandSides[rowIdx],
                                   *upperBound))
      tryUpdateBest(rowIdx,
                    _simplexTableau._rightHandSides[rowIdx] - *upperBound);
  }

  return bestRowIdx;
}

template <typename T, typename SimplexTraitsT>
std::optional<int>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>::chooseEnteringColumnIdx(
    const int pivotRowIdx, const std::vector<T> &pivotRow,
    const bool isPivotRowUnderLowerBound) {
  std::optional<int> mostRestrictiveColumnIdx;
  std::optional<T> mostRestrictiveColumnBound;

  auto &isColumnAtLowerBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset;

  const auto currentRestrictionBound =
      [&](const int colIdx) -> std::optional<T> {
    if (SimplexTraitsT::isEligibleForPivot(pivotRow[colIdx])) {
      if (isPivotRowUnderLowerBound) {
        if (isColumnAtLowerBoundBitset[colIdx]) {
          if (SimplexTraitsT::less(pivotRow[colIdx], 0.0))
            return _simplexTableau._reducedCosts[colIdx] / (-pivotRow[colIdx]);
        } else {
          if (SimplexTraitsT::greater(pivotRow[colIdx], 0.0))
            return (-_simplexTableau._reducedCosts[colIdx]) / pivotRow[colIdx];
        }
      } else {
        if (isColumnAtLowerBoundBitset[colIdx]) {
          if (SimplexTraitsT::greater(pivotRow[colIdx], 0.0))
            return _simplexTableau._reducedCosts[colIdx] / pivotRow[colIdx];
        } else {
          if (SimplexTraitsT::less(pivotRow[colIdx], 0.0))
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

template class RevisedDualSimplexPFIBounds<double>;
template class RevisedDualSimplexPFIBounds<long double>;