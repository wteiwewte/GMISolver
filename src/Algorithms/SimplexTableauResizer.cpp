#include "src/Algorithms/SimplexTableauResizer.h"

#include "src/Algorithms/SimplexTableau.h"

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
SimplexTableauResizer<T, SimplexTraitsT>::SimplexTableauResizer(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager) {
}

template <typename T, typename SimplexTraitsT>
void SimplexTableauResizer<
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
void SimplexTableauResizer<
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
void SimplexTableauResizer<T, SimplexTraitsT>::removeRows(
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

template class SimplexTableauResizer<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableauResizer<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class SimplexTableauResizer<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableauResizer<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;