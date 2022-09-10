#include "src/Algorithms/RevisedPrimalSimplexPFI.h"

#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/SpdlogHeader.h"


template <typename T, typename ComparisonTraitsT>
RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::RevisedPrimalSimplexPFI(
    SimplexTableau<T> &simplexTableau)
    : _simplexTableau(simplexTableau) {}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::runPhaseOne() {
  SPDLOG_INFO("BASIS SIZE {}", _simplexTableau._rowInfos.size());
  run();
  removeArtificialVariablesFromBasis();
  removeArtificialVariablesFromProgram();

  if (ComparisonTraitsT::greater(_simplexTableau._objectiveValue, 0.0)) {
    SPDLOG_INFO("Program with artificial variable has optimum greater than 0 - "
                "initial program is infeasible");
    return;
  }

  setInitialObjective();
  run();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::run() {
  int iterCount = 1;
  while (true) {
    SPDLOG_INFO("ITERATION {}", iterCount++);
    const bool iterResult = runOneIteration();
    if (iterResult)
      break;

    _simplexTableau.calculateCurrentObjectiveValue();
    _simplexTableau.calculateSolution();
    SPDLOG_INFO("{}\n", _simplexTableau.toStringShort());

    //    if (iterCount < 10)
    //      SPDLOG_INFO("{}\n", _simplexTableau.toString());

    //    reinversion();
    //    _simplexTableau.calculateCurrentObjectiveValue();
    //    _simplexTableau.calculateSolution();

    //    if (iterCount < 10)
    //      SPDLOG_INFO("{}\n", _simplexTableau.toString());
  }
  //  SPDLOG_INFO("{}\n", _simplexTableau.toString());
  SPDLOG_INFO("{}\n", _simplexTableau.toStringShort());
  SPDLOG_INFO("SIMPLEX ENDED");
}
template <typename T, typename ComparisonTraitsT>
bool RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::runOneIteration() {
  std::optional<int> enteringColumnIdx;
  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size();
       ++columnIdx) {
    if (_simplexTableau._variableInfos[columnIdx]._isArtificial)
      continue;

    if (!_simplexTableau._simplexBasisData
             ._isBasicColumnIndexBitset[columnIdx] &&
        ComparisonTraitsT::less(_simplexTableau._reducedCosts[columnIdx],
                                0.0)) {
      enteringColumnIdx = columnIdx;
      SPDLOG_DEBUG("ENTERING COLUMN IDX {} REDUCED COST {}", columnIdx,
                   _simplexTableau._reducedCosts[columnIdx]);

      const std::vector<T> enteringColumn =
          computeEnteringColumn(*enteringColumnIdx);
      const std::optional<int> rowIdx = chooseRowIdx(enteringColumn);
      if (!rowIdx.has_value()) {
        _simplexTableau._result = LPOptimizationResult::UNBOUNDED;
        return true;
      }
      const std::vector<T> pivotRow = computePivotRow(*rowIdx);

      pivot(*rowIdx, *enteringColumnIdx, enteringColumn, pivotRow);
      return false;
    }
  }

  _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
  return true;
}

template <typename T, typename ComparisonTraitsT>
std::vector<T>
RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::computeEnteringColumn(
    const int enteringColumnIdx) {
  std::vector<T> result(_simplexTableau._rowInfos.size());

  // TODO - maybe add kahan summation algo, maybe opt order
  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i) {
    for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k)
      result[i] += _simplexTableau._basisMatrixInverse[i][k] *
                   _simplexTableau._constraintMatrix[k][enteringColumnIdx];
  }

  return result;
}

template <typename T, typename ComparisonTraitsT>
std::vector<T> RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::computePivotRow(
    const int rowIdx) {
  std::vector<T> result(_simplexTableau._variableInfos.size());

  // TODO - maybe add kahan summation algo, maybe opt order
  const int columnBasicIdx =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
    if (j == columnBasicIdx)
      result[j] = 1.0;
    else if (j != columnBasicIdx &&
             _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
      result[j] = 0.0;
    else
      for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k)
        result[j] += _simplexTableau._basisMatrixInverse[rowIdx][k] *
                     _simplexTableau._constraintMatrix[k][j];
  }

  return result;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::pivot(
    const int rowIdx, const int enteringColumnIdx,
    const std::vector<T> &enteringColumn, const std::vector<T> &pivotRow) {
  const PivotData<T> pivotData{rowIdx, enteringColumnIdx,
                               1.0 / enteringColumn[rowIdx]};
  SPDLOG_INFO("PIVOT COEFF {}, INV {}", enteringColumn[rowIdx],
              1.0 / enteringColumn[rowIdx]);
  SPDLOG_INFO("RHS {}", _simplexTableau._rightHandSides[rowIdx]);
  updateReducedCosts(pivotData, pivotRow);
  updateInverseMatrixWithRHS(pivotData, enteringColumn);
  _simplexTableau.updateBasisData(pivotData);
}

template <typename T, typename ComparisonTraitsT>
std::optional<int> RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::chooseRowIdx(
    const std::vector<T> &enteringColumn) {
  std::optional<int> leavingRowIdx;

  // TODO - add choosing most negative

  for (int rowIdx = 0; rowIdx < _simplexTableau.getRowInfos().size(); ++rowIdx)
    if (ComparisonTraitsT::greater(enteringColumn[rowIdx], 0.0))
      if (!leavingRowIdx.has_value() ||
          ComparisonTraitsT::less(
              _simplexTableau._rightHandSides[rowIdx] *
                  enteringColumn[*leavingRowIdx],
              _simplexTableau._rightHandSides[*leavingRowIdx] *
                  enteringColumn[rowIdx]))
        leavingRowIdx = rowIdx;

  return leavingRowIdx;
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::updateReducedCosts(
    const PivotData<T> &pivotData, const std::vector<T> &pivotRow) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;
  const auto commonCoeffReducedCost =
      _simplexTableau._reducedCosts[enteringColumnIdx] * pivotingTermInverse;
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    _simplexTableau._reducedCosts[j] -= commonCoeffReducedCost * pivotRow[j];
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::updateInverseMatrixWithRHS(
    const PivotData<T> &pivotData, const std::vector<T> &enteringColumn) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;

  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i) {
    if (i == leavingRowIdx)
      continue;
    // TODO - opt if coeff is zero
    const auto commonCoeff = enteringColumn[i] * pivotingTermInverse;

    for (int j = 0; j < _simplexTableau._basisMatrixInverse.size(); ++j)
      _simplexTableau._basisMatrixInverse[i][j] -=
          commonCoeff * _simplexTableau._basisMatrixInverse[leavingRowIdx][j];

    _simplexTableau._rightHandSides[i] -=
        commonCoeff * _simplexTableau._rightHandSides[leavingRowIdx];
  }

  for (int j = 0; j < _simplexTableau._basisMatrixInverse.size(); ++j)
    _simplexTableau._basisMatrixInverse[leavingRowIdx][j] *=
        pivotingTermInverse;

  _simplexTableau._rightHandSides[leavingRowIdx] *= pivotingTermInverse;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<
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
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    _simplexTableau._constraintMatrix[rowIdx].resize(*firstArtificialIdx);
  }
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<
    T, ComparisonTraitsT>::removeArtificialVariablesFromBasis() {
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicVarIdx =
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    SPDLOG_INFO("FOUND BASIC ARTIFICIAL COLUMN IDX {}", basicVarIdx);

    std::optional<int> nonZeroEntryColumnIndex;
    const std::vector<T> pivotRow = computePivotRow(rowIdx);

    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j) {
      if (!_simplexTableau._variableInfos[j]._isArtificial ||
          _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        continue;

      if (!ComparisonTraitsT::equal(pivotRow[j], 0.0)) {
        nonZeroEntryColumnIndex = j;
        break;
      }
    }

    if (!nonZeroEntryColumnIndex.has_value()) {
      SPDLOG_WARN("Redundant constraints in lp formulation!");
      //      removeRow(rowIdx);
    } else {
      const std::vector<T> enteringColumn =
          computeEnteringColumn(*nonZeroEntryColumnIndex);
      pivot(rowIdx, *nonZeroEntryColumnIndex, enteringColumn, pivotRow);
    }
  }
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::removeRow(
    const int rowIdx) {
  // rethink it
  //  auto& [rowToBasisColumnIdxMap, isBasicColumnIndexBitset] =
  //  _simplexTableau._simplexBasisData;
  //  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[rowIdx]] = false;
  //  for (int i = rowIdx; i < _simplexTableau._rowInfos.size() - 1; ++i)
  //    rowToBasisColumnIdxMap[i] = rowToBasisColumnIdxMap[i + 1];
  //  rowToBasisColumnIdxMap.pop_back();
  //
  //  _simplexTableau._constraintMatrix.erase(_simplexTableau._constraintMatrix.begin()
  //  + rowIdx);
  //  _simplexTableau._rowInfos.erase(_simplexTableau._rowInfos.begin() +
  //  rowIdx);
  //  _simplexTableau._rightHandSides.erase(_simplexTableau._rightHandSides.begin()
  //  + rowIdx); for (int i = 0; i < _simplexTableau._basisMatrixInverse.size();
  //  ++i)
  //    _simplexTableau._basisMatrixInverse[i].erase(_simplexTableau._basisMatrixInverse[i].begin()
  //    + rowIdx);
  //  _simplexTableau._basisMatrixInverse.erase(_simplexTableau._basisMatrixInverse.begin()
  //  + rowIdx);
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::setInitialObjective() {
  _simplexTableau._objectiveRow =
      _simplexTableau._initialProgram.getObjective();
  calculateDual();
  _simplexTableau.calculateReducedCostsBasedOnDual();
  _simplexTableau.calculateCurrentObjectiveValue();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::calculateDual() {
  _simplexTableau._y.resize(_simplexTableau._basisMatrixInverse.size());
  for (int colIdx = 0; colIdx < _simplexTableau._y.size(); ++colIdx) {
    T sum{};
    for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k)
      sum += _simplexTableau._objectiveRow[_simplexTableau._simplexBasisData
                                               ._rowToBasisColumnIdxMap[k]] *
             _simplexTableau._basisMatrixInverse[k][colIdx];

    _simplexTableau._y[colIdx] = sum;
  }
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::reinversion() {
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
  std::vector<T> newRHS = _simplexTableau._initialProgram.getRightHandSides();
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    newBasisMatrixInverse[rowIdx].resize(basisSize);
    newBasisMatrixInverse[rowIdx][rowIdx] = 1.0;
  }

  std::vector<bool> isUnusedColumn(basisSize, true);

  const auto findPivotColumn = [&](const int rowIdx) -> std::optional<int> {
    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx] &&
          !ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
        return colIdx;

    return std::nullopt;
  };

  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    const auto pivotColumnIdx = findPivotColumn(rowIdx);
    if (!pivotColumnIdx.has_value()) {
      SPDLOG_WARN("Basis matrix reinversion failed for row {}!", rowIdx);
      return;
    }

    const T pivotingTermInverse{1.0 / basisColumns[*pivotColumnIdx][rowIdx]};
    //    SPDLOG_DEBUG("REINVERSION PIVOT {}", pivotingTermInverse);
    for (int j = 0; j < basisSize; ++j) {
      if (j == rowIdx)
        continue; // !!

      if (ComparisonTraitsT::equal(basisColumns[*pivotColumnIdx][j], 0.0))
        continue;

      const auto commonCoeff =
          pivotingTermInverse * basisColumns[*pivotColumnIdx][j];
      for (int k = 0; k < basisSize; ++k)
        newBasisMatrixInverse[j][k] -=
            commonCoeff * newBasisMatrixInverse[rowIdx][k];

      newRHS[j] -= commonCoeff * newRHS[rowIdx];
    }

    for (int k = 0; k < basisSize; ++k)
      newBasisMatrixInverse[rowIdx][k] *= pivotingTermInverse;

    newRHS[rowIdx] *= pivotingTermInverse;

    isUnusedColumn[*pivotColumnIdx] = false;
    _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
        columnIndexesMapping[*pivotColumnIdx];

    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx]) {
        if (ComparisonTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
          continue;

        for (int k = 0; k < basisSize; ++k) {
          if (k == rowIdx)
            continue;

          basisColumns[colIdx][k] -= pivotingTermInverse *
                                     basisColumns[*pivotColumnIdx][k] *
                                     basisColumns[colIdx][rowIdx];
        }

        basisColumns[colIdx][rowIdx] *= pivotingTermInverse;
      }
  }

  _simplexTableau._basisMatrixInverse = newBasisMatrixInverse;
  _simplexTableau._rightHandSides = newRHS;
  SPDLOG_INFO("REINVERSION SUCCESS");
}
// template <typename T, typename ComparisonTraitsT>
// void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::calculateRHS() {
//   for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
//   {
//
//   }
// }

template class RevisedPrimalSimplexPFI<double>;