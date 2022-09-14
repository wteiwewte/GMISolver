#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"

#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/SpdlogHeader.h"

template <typename T, typename ComparisonTraitsT>
RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::RevisedDualSimplexPFIBounds(
    SimplexTableau<T> &simplexTableau,
    const DualSimplexRowPivotRule dualSimplexRowPivotRule,
    const int32_t objValueLoggingFrequency, const int32_t reinversionFrequency)
    : _simplexTableau(simplexTableau),
      _dualSimplexRowPivotRule(dualSimplexRowPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency) {
  initRHS();
  calculateCurrentObjectiveValue();
  calculateSolution();
}

template <typename T, typename ComparisonTraitsT>
void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::run() {
  SPDLOG_INFO("BASIS SIZE {} COLUMN PIVOT RULE {}",
              _simplexTableau._rowInfos.size(),
              dualSimplexRowPivotRuleToStr(_dualSimplexRowPivotRule));
  SPDLOG_TRACE("{}\n", _simplexTableau.toString());

  [[maybe_unused]] int iterCount = 1;
  while (true) {
    const bool iterResult = runOneIteration();
    if (iterResult)
      break;

    calculateCurrentObjectiveValue();
    calculateSolution();

    ++iterCount;
    if (_objValueLoggingFrequency &&
        (iterCount % _objValueLoggingFrequency == 0)) {
      SPDLOG_INFO("ITERATION {}", iterCount);
      SPDLOG_INFO("{}\n", _simplexTableau.toStringObjectiveValue());
    }
    SPDLOG_DEBUG("{}\n", _simplexTableau.toStringObjectiveValue());
  }
  SPDLOG_INFO("DUAL SIMPLEX ENDED, ITERATION COUNT {}", iterCount);
}

template <typename T, typename ComparisonTraitsT>
bool RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::runOneIteration() {
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
  const bool isPivotRowUnderLowerBound = ComparisonTraitsT::less(
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

  pivot(pivotRow, *pivotRowIdx, enteringColumn, *enteringColumnIdx,
        isPivotRowUnderLowerBound);
  return false;
}

template <typename T, typename ComparisonTraitsT>
std::optional<int>
RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::chooseRow() {
  switch (_dualSimplexRowPivotRule) {
  case DualSimplexRowPivotRule::FIRST_ELIGIBLE:
    return chooseRowFirstEligible();
  case DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION:
    return chooseRowBiggestViolation();
  }
}

template <typename T, typename ComparisonTraitsT>
std::optional<int>
RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::chooseRowFirstEligible() {
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicColumnIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    const auto lowerBound =
        _simplexTableau._variableLowerBounds[basicColumnIdx];
    if (lowerBound.has_value() &&
        ComparisonTraitsT::less(_simplexTableau._rightHandSides[rowIdx],
                                *lowerBound))
      return rowIdx;

    const auto upperBound =
        _simplexTableau._variableUpperBounds[basicColumnIdx];
    if (upperBound.has_value() &&
        ComparisonTraitsT::greater(_simplexTableau._rightHandSides[rowIdx],
                                   *upperBound))
      return rowIdx;
  }
  return std::nullopt;
}

template <typename T, typename ComparisonTraitsT>
std::optional<int>
RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::chooseRowBiggestViolation() {
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
        ComparisonTraitsT::less(_simplexTableau._rightHandSides[rowIdx],
                                *lowerBound))
      tryUpdateBest(rowIdx,
                    *lowerBound - _simplexTableau._rightHandSides[rowIdx]);

    const auto upperBound =
        _simplexTableau._variableUpperBounds[basicColumnIdx];
    if (upperBound.has_value() &&
        ComparisonTraitsT::greater(_simplexTableau._rightHandSides[rowIdx],
                                   *upperBound))
      tryUpdateBest(rowIdx,
                    _simplexTableau._rightHandSides[rowIdx] - *upperBound);
  }

  return bestRowIdx;
}

template <typename T, typename ComparisonTraitsT>
std::optional<int>
RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::chooseEnteringColumnIdx(
    const int pivotRowIdx, const std::vector<T> &pivotRow,
    const bool isPivotRowUnderLowerBound) {
  std::optional<int> mostRestrictiveColumnIdx;
  std::optional<T> mostRestrictiveColumnBound;

  auto &isColumnAtLowerBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset;

  const auto currentRestrictionBound =
      [&](const int colIdx) -> std::optional<T> {
    if (ComparisonTraitsT::isEligibleForPivot(pivotRow[colIdx])) {
      if (isPivotRowUnderLowerBound) {
        if (isColumnAtLowerBoundBitset[colIdx]) {
          if (ComparisonTraitsT::less(pivotRow[colIdx], 0.0))
            return _simplexTableau._reducedCosts[colIdx] / (-pivotRow[colIdx]);
        } else {
          if (ComparisonTraitsT::greater(pivotRow[colIdx], 0.0))
            return (-_simplexTableau._reducedCosts[colIdx]) / pivotRow[colIdx];
        }
      } else {
        if (isColumnAtLowerBoundBitset[colIdx]) {
          if (ComparisonTraitsT::greater(pivotRow[colIdx], 0.0))
            return _simplexTableau._reducedCosts[colIdx] / pivotRow[colIdx];
        } else {
          if (ComparisonTraitsT::less(pivotRow[colIdx], 0.0))
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

template <typename T, typename ComparisonTraitsT>
void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::pivot(
    const std::vector<T> &pivotRow, const int pivotRowIdx,
    const std::vector<T> &enteringColumn, const int enteringColumnIdx,
    const bool isPivotRowUnderLowerBound) {
  auto &isColumnAtLowerBoundBitset =
      _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset;

  SPDLOG_DEBUG("PIVOT - ENTERING COLUMN IDX {}, ROW IDX {}", enteringColumnIdx,
               pivotRowIdx);
  const auto leavingBasicColumnIdx =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[pivotRowIdx];
  const auto leavingColumn =
      _simplexTableau.computeTableauColumn(leavingBasicColumnIdx);

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    _simplexTableau._rightHandSides[rowIdx] +=
        (isColumnAtLowerBoundBitset[enteringColumnIdx]
             ? *_simplexTableau._variableLowerBounds[enteringColumnIdx]
             : *_simplexTableau._variableUpperBounds[enteringColumnIdx]) *
        enteringColumn[rowIdx];

    _simplexTableau._rightHandSides[rowIdx] -=
        (isPivotRowUnderLowerBound
             ? *_simplexTableau._variableLowerBounds[leavingBasicColumnIdx]
             : *_simplexTableau._variableUpperBounds[leavingBasicColumnIdx]) *
        leavingColumn[rowIdx];
  }

  _simplexTableau.pivotImplicitBounds(pivotRowIdx, enteringColumnIdx,
                                      enteringColumn, pivotRow,
                                      isPivotRowUnderLowerBound);
}

//
// template <typename T, typename ComparisonTraitsT>
// std::optional<PivotRowData<T>>
// RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::chooseRowIdx(
//    const int enteringColumnIdx, const std::vector<T> &enteringColumn) {
//  std::optional<int> rowIdxMostRestrictiveLowerBound;
//  std::optional<T> mostRestrictiveLowerBound;
//
//  std::optional<int> rowIdxMostRestrictiveUpperBound;
//  std::optional<T> mostRestrictiveUpperBound;
//
//  for (int rowIdx = 0; rowIdx < _simplexTableau.getRowInfos().size();
//       ++rowIdx) {
//    if (ComparisonTraitsT::equal(enteringColumn[rowIdx], 0.0))
//      continue;
//
//    const auto basicColumnIdx =
//        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
//    const auto lowerBound =
//        _simplexTableau._variableLowerBounds[basicColumnIdx];
//    const auto upperBound =
//        _simplexTableau._variableUpperBounds[basicColumnIdx];
//
//    if (_simplexTableau._simplexBasisData
//            ._isColumnAtLowerBoundBitset[enteringColumnIdx]) {
//      if (enteringColumn[rowIdx] > 0.0) {
//        if (!lowerBound.has_value())
//          continue;
//
//        const auto currentRestrictionValue =
//            (_simplexTableau._rightHandSides[rowIdx] - *lowerBound) /
//            enteringColumn[rowIdx];
//        if (!rowIdxMostRestrictiveLowerBound.has_value() ||
//            (currentRestrictionValue < *mostRestrictiveLowerBound)) {
//          rowIdxMostRestrictiveLowerBound = rowIdx;
//          mostRestrictiveLowerBound = currentRestrictionValue;
//        }
//      } else {
//        if (!upperBound.has_value())
//          continue;
//
//        const auto currentRestrictionValue =
//            (*upperBound - _simplexTableau._rightHandSides[rowIdx]) /
//            (-enteringColumn[rowIdx]);
//        if (!rowIdxMostRestrictiveUpperBound.has_value() ||
//            (currentRestrictionValue < *mostRestrictiveUpperBound)) {
//          rowIdxMostRestrictiveUpperBound = rowIdx;
//          mostRestrictiveUpperBound = currentRestrictionValue;
//        }
//      }
//    } else {
//      if (enteringColumn[rowIdx] > 0.0) {
//        if (!upperBound.has_value())
//          continue;
//
//        const auto currentRestrictionValue =
//            (*upperBound - _simplexTableau._rightHandSides[rowIdx]) /
//            enteringColumn[rowIdx];
//        if (!rowIdxMostRestrictiveUpperBound.has_value() ||
//            (currentRestrictionValue < *mostRestrictiveUpperBound)) {
//          rowIdxMostRestrictiveUpperBound = rowIdx;
//          mostRestrictiveUpperBound = currentRestrictionValue;
//        }
//      } else {
//        if (!lowerBound.has_value())
//          continue;
//
//        const auto currentRestrictionValue =
//            (_simplexTableau._rightHandSides[rowIdx] - *lowerBound) /
//            (-enteringColumn[rowIdx]);
//        if (!rowIdxMostRestrictiveLowerBound.has_value() ||
//            (currentRestrictionValue < *mostRestrictiveLowerBound)) {
//          rowIdxMostRestrictiveLowerBound = rowIdx;
//          mostRestrictiveLowerBound = currentRestrictionValue;
//        }
//      }
//    }
//  }
//
//  PivotRowData<T> pivotRowData;
//  const auto enteringColumnLowerBound =
//      _simplexTableau._variableLowerBounds[enteringColumnIdx];
//  const auto enteringColumnUpperBound =
//      _simplexTableau._variableUpperBounds[enteringColumnIdx];
//  if (enteringColumnLowerBound.has_value() &&
//      enteringColumnUpperBound.has_value()) {
//    pivotRowData._minRatio =
//        *enteringColumnUpperBound - *enteringColumnLowerBound;
//    pivotRowData._noBasisChangeNeeded = true;
//  }
//
//  if (rowIdxMostRestrictiveLowerBound.has_value()) {
//    if (!pivotRowData._minRatio.has_value() ||
//        (*mostRestrictiveLowerBound < pivotRowData._minRatio)) {
//      pivotRowData._minRatio = *mostRestrictiveLowerBound;
//      pivotRowData._pivotRowIdx = *rowIdxMostRestrictiveLowerBound;
//      pivotRowData._departingIdxBecomesLowerBound = true;
//      pivotRowData._noBasisChangeNeeded = false;
//    }
//  }
//
//  if (rowIdxMostRestrictiveUpperBound.has_value()) {
//    if (!pivotRowData._minRatio.has_value() ||
//        (*mostRestrictiveUpperBound < pivotRowData._minRatio)) {
//      pivotRowData._minRatio = *mostRestrictiveUpperBound;
//      pivotRowData._pivotRowIdx = *rowIdxMostRestrictiveUpperBound;
//      pivotRowData._departingIdxBecomesLowerBound = false;
//      pivotRowData._noBasisChangeNeeded = false;
//    }
//  }
//
//  return pivotRowData._minRatio.has_value() ? std::optional{pivotRowData}
//                                            : std::nullopt;
//}
//
//
// template <typename T, typename ComparisonTraitsT>
// void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::calculateDual() {
//  _simplexTableau._y.resize(_simplexTableau._basisMatrixInverse.size());
//  for (int colIdx = 0; colIdx < _simplexTableau._y.size(); ++colIdx) {
//    T sum{};
//
//    for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k) {
//      SPDLOG_INFO("{} {} {}", k,
//                   _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[k],
//                   _simplexTableau._objectiveRow.size());
//      sum += _simplexTableau._objectiveRow[_simplexTableau._simplexBasisData
//                                               ._rowToBasisColumnIdxMap[k]] *
//             _simplexTableau._basisMatrixInverse[k][colIdx];
//    }
//
//    _simplexTableau._y[colIdx] = sum;
//  }
//}
//
template <typename T, typename ComparisonTraitsT>
void RevisedDualSimplexPFIBounds<
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
//
//// template <typename T, typename ComparisonTraitsT>
//// void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::calculateSolution()
//// {
////   _simplexTableau._x.resize(_simplexTableau._variableInfos.size());
////   std::fill(_simplexTableau._x.begin(), _simplexTableau._x.end(), T{});
////   for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
////
///_simplexTableau._x[_simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx]]
////     =
////         _simplexTableau._rightHandSides[rowIdx];
////
////   for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
////     if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
////       _simplexTableau._x[j]  =
////       (_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[j] ?
//// *_simplexTableau._variableLowerBounds[j] :
///*_simplexTableau._variableUpperBounds[j]); / }
template <typename T, typename ComparisonTraitsT>
void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::calculateSolution() {
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
//
template <typename T, typename ComparisonTraitsT>
void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::initRHS() {
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
// template <typename T, typename ComparisonTraitsT>
// void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::calculateRHS() {
//   std::vector<T> tempRHS = _simplexTableau._initialRightHandSides;
//   for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
//     //    const auto basicColumnIdx = _simplexTableau._simplexB
//     //    asisData._rowToBasisColumnIdxMap[rowIdx];
//     for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
//       if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
//         tempRHS[rowIdx] -=
//             (_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[j]
//                  ? *_simplexTableau._variableLowerBounds[j]
//                  : *_simplexTableau._variableUpperBounds[j]) *
//             _simplexTableau._constraintMatrix[rowIdx][j];
//   }
//
//   std::fill(_simplexTableau._rightHandSides.begin(),
//             _simplexTableau._rightHandSides.end(), T{0.0});
//   // TODO - summation
//   for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i)
//     for (int j = 0; j < _simplexTableau._rowInfos.size(); ++j)
//       _simplexTableau._rightHandSides[i] +=
//           _simplexTableau._basisMatrixInverse[i][j] * tempRHS[j];
// }
//
// template <typename T, typename ComparisonTraitsT>
// void RevisedDualSimplexPFIBounds<T, ComparisonTraitsT>::reinversion() {
//   // TODO - opt it by storing columns in vectors
//   const auto basisSize = _simplexTableau._rowInfos.size();
//   std::vector<int> columnIndexesMapping(basisSize);
//   std::vector<std::vector<T>> basisColumns(basisSize);
//   for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
//     const auto basicColumnIdx =
//         _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
//     columnIndexesMapping[rowIdx] = basicColumnIdx;
//     basisColumns[rowIdx].resize(basisSize);
//     for (int j = 0; j < basisSize; ++j)
//       basisColumns[rowIdx][j] =
//           _simplexTableau._constraintMatrix[j][basicColumnIdx];
//   }
//
//   std::vector<std::vector<T>> newBasisMatrixInverse(basisSize);
//   std::vector<T> newRHS =
//   _simplexTableau._initialProgram.getRightHandSides(); for (int rowIdx = 0;
//   rowIdx < basisSize; ++rowIdx) {
//     newBasisMatrixInverse[rowIdx].resize(basisSize);
//     newBasisMatrixInverse[rowIdx][rowIdx] = 1.0;
//   }
//
//   std::vector<bool> isUnusedColumn(basisSize, true);
//
//   const auto findPivotColumn = [&](const int rowIdx) -> std::optional<int> {
//     for (int colIdx = 0; colIdx < basisSize; ++colIdx)
//       if (isUnusedColumn[colIdx] &&
//           !ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
//         return colIdx;
//
//     return std::nullopt;
//   };
//
//   for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
//     const auto pivotColumnIdx = findPivotColumn(rowIdx);
//     if (!pivotColumnIdx.has_value()) {
//       SPDLOG_WARN("Basis matrix reinversion failed for row {}!", rowIdx);
//       return;
//     }
//
//     const T pivotingTermInverse{1.0 / basisColumns[*pivotColumnIdx][rowIdx]};
//     //    SPDLOG_DEBUG("REINVERSION PIVOT {}", pivotingTermInverse);
//     for (int j = 0; j < basisSize; ++j) {
//       if (j == rowIdx)
//         continue; // !!
//
//       if (ComparisonTraitsT::equal(basisColumns[*pivotColumnIdx][j], 0.0))
//         continue;
//
//       const auto commonCoeff =
//           pivotingTermInverse * basisColumns[*pivotColumnIdx][j];
//       for (int k = 0; k < basisSize; ++k)
//         newBasisMatrixInverse[j][k] -=
//             commonCoeff * newBasisMatrixInverse[rowIdx][k];
//
//       newRHS[j] -= commonCoeff * newRHS[rowIdx];
//     }
//
//     for (int k = 0; k < basisSize; ++k)
//       newBasisMatrixInverse[rowIdx][k] *= pivotingTermInverse;
//
//     newRHS[rowIdx] *= pivotingTermInverse;
//
//     isUnusedColumn[*pivotColumnIdx] = false;
//     _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
//         columnIndexesMapping[*pivotColumnIdx];
//
//     for (int colIdx = 0; colIdx < basisSize; ++colIdx)
//       if (isUnusedColumn[colIdx]) {
//         if (ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
//           continue;
//
//         for (int k = 0; k < basisSize; ++k) {
//           if (k == rowIdx)
//             continue;
//
//           basisColumns[colIdx][k] -= pivotingTermInverse *
//                                      basisColumns[*pivotColumnIdx][k] *
//                                      basisColumns[colIdx][rowIdx];
//         }
//
//         basisColumns[colIdx][rowIdx] *= pivotingTermInverse;
//       }
//   }
//
//   _simplexTableau._basisMatrixInverse = newBasisMatrixInverse;
//   //  _simplexTableau._rightHandSides = newRHS;
//   calculateRHS();
//   SPDLOG_INFO("REINVERSION SUCCESS");
// }

template class RevisedDualSimplexPFIBounds<double>;