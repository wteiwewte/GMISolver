#include "src/Algorithms/SimplexTableauResizer.h"

#include "src/Algorithms/ReinversionManager.h"
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

void removeElements(boost::dynamic_bitset<> &bitset,
                    const std::vector<bool> &shouldBeRemoved) {
  auto copyBitset = bitset;
  int currentNewIndex = 0;
  for (int i = 0; i < bitset.size(); ++i) {
    if (!shouldBeRemoved[i])
      copyBitset[currentNewIndex++] = bitset[i];
  }
  copyBitset.resize(currentNewIndex);
  bitset = copyBitset;
}
} // namespace

template <typename T, typename SimplexTraitsT>
SimplexTableauResizer<T, SimplexTraitsT>::SimplexTableauResizer(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager) {
}

template <typename T, typename SimplexTraitsT>
bool SimplexTableauResizer<
    T, SimplexTraitsT>::removeArtificialVariablesFromProgram() {
  const std::vector<bool> shouldRowBeRemoved =
      moveArtificialVariablesOutOfBasis();
  removeRows(shouldRowBeRemoved);

  std::vector<bool> shouldVarBeRemoved(_simplexTableau._variableInfos.size(),
                                       false);
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
    if (_simplexTableau._variableInfos[j]._isArtificial)
      shouldVarBeRemoved[j] = true;
  }

  removeVariables(shouldVarBeRemoved);
  _simplexTableau.initMatrixRepresentations();
  return _reinversionManager.reinverse();
}

template <typename T, typename SimplexTraitsT>
std::vector<bool>
SimplexTableauResizer<T, SimplexTraitsT>::moveArtificialVariablesOutOfBasis() {
  std::vector<bool> shouldRowBeRemoved(_simplexTableau._rowInfos.size(), false);

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    SPDLOG_DEBUG("FOUND BASIC ARTIFICIAL ROW IDX {} COLUMN IDX {}, RHS {}",
                 rowIdx, basicVarIdx, _simplexTableau._rightHandSides[rowIdx]);

    std::optional<int> nonZeroEntryColumnIndex;
    const auto pivotRowSharedPtr =
        _simplexTableau.computeTableauRowGeneric(rowIdx);
    const auto &pivotRow = *pivotRowSharedPtr;

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

  return shouldRowBeRemoved;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableauResizer<T, SimplexTraitsT>::removeRows(
    const std::vector<bool> &shouldRowBeRemoved) {
  std::vector<int> rowsToBeRemoved;
  for (int rowIdx = 0; rowIdx < shouldRowBeRemoved.size(); ++rowIdx)
    if (shouldRowBeRemoved[rowIdx])
      rowsToBeRemoved.push_back(rowIdx);

  if (rowsToBeRemoved.size() == 0)
    return;

  SPDLOG_INFO("{} CONSTRAINTS [{}] TO BE REMOVED", rowsToBeRemoved.size(),
              fmt::join(rowsToBeRemoved, ", "));
  removeElements(_simplexTableau._rowInfos, shouldRowBeRemoved);
  removeElements(_simplexTableau._rightHandSides, shouldRowBeRemoved);
  removeElements(_simplexTableau._simplexBasisData._rowToBasisColumnIdxMap,
                 shouldRowBeRemoved);

  switch (_simplexTableau._simplexTableauType) {
  case SimplexTableauType::FULL: {
    removeElements(_simplexTableau._fullTableau, shouldRowBeRemoved);
    break;
  }
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE: {
    removeElements(_simplexTableau._constraintMatrix, shouldRowBeRemoved);
    removeElements(_simplexTableau._initialRightHandSides, shouldRowBeRemoved);
    for (int i = 0; i < _simplexTableau._basisMatrixInverse.size(); ++i) {
      removeElements(_simplexTableau._basisMatrixInverse[i],
                     shouldRowBeRemoved);
    }
    removeElements(_simplexTableau._basisMatrixInverse, shouldRowBeRemoved);
    break;
  }
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
    removeElements(_simplexTableau._constraintMatrix, shouldRowBeRemoved);
    removeElements(_simplexTableau._initialRightHandSides, shouldRowBeRemoved);
    break;
  }
  default:
    break;
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableauResizer<T, SimplexTraitsT>::removeVariables(
    const std::vector<bool> &shouldVarBeRemoved) {
  removeElements(_simplexTableau._variableInfos, shouldVarBeRemoved);
  removeElements(_simplexTableau._isVariableFreeBitset, shouldVarBeRemoved);
  removeElements(_simplexTableau._variableLowerBounds, shouldVarBeRemoved);
  removeElements(_simplexTableau._variableUpperBounds, shouldVarBeRemoved);

  removeElements(_simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset,
                 shouldVarBeRemoved);
  removeElements(_simplexTableau._simplexBasisData._isColumnAtUpperBoundBitset,
                 shouldVarBeRemoved);
  removeElements(_simplexTableau._simplexBasisData._isBasicColumnIndexBitset,
                 shouldVarBeRemoved);

  removeElements(_simplexTableau._reducedCosts, shouldVarBeRemoved);
  removeElements(_simplexTableau._objectiveRow, shouldVarBeRemoved);
  removeElements(_simplexTableau._x, shouldVarBeRemoved);
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    if (_simplexTableau._simplexTableauType == SimplexTableauType::FULL) {
      removeElements(_simplexTableau._fullTableau[rowIdx], shouldVarBeRemoved);
    } else {
      removeElements(_simplexTableau._constraintMatrix[rowIdx],
                     shouldVarBeRemoved);
    }
  }
}

template class SimplexTableauResizer<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableauResizer<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class SimplexTableauResizer<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableauResizer<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;