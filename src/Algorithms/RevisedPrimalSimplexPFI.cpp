#include "src/Algorithms/RevisedPrimalSimplexPFI.h"

#include "src/Algorithms/SimplexTableau.h"

template <typename T, typename ComparisonTraitsT>
RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::RevisedPrimalSimplexPFI(
    SimplexTableau<T> &simplexTableau) : _simplexTableau(simplexTableau) {}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::runPhaseOne() {
  run();
  removeArtificialVariablesFromBasis();
  removeArtificialVariablesFromProgram();

  if (ComparisonTraitsT::greater(_simplexTableau._objectiveValue, 0.0))
  {
    spdlog::info("Program with artificial variable has optimum greater than 0 - "
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
    spdlog::info("ITERATION {}", iterCount++);
    const bool iterResult = runOneIteration();
    if (iterResult)
      break;

    _simplexTableau.calculateCurrentObjectiveValue();
    _simplexTableau.calculateSolution();
//    spdlog::info("{}\n", _simplexTableau.toString());
//    spdlog::info("{}\n", _simplexTableau.toStringShort());
  }
  spdlog::info("SIMPLEX ENDED");
}
template <typename T, typename ComparisonTraitsT>
bool RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::runOneIteration() {
  std::optional<int> enteringColumnIdx;
  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size(); ++columnIdx) {
    if(_simplexTableau._variableInfos[columnIdx]._isArtificial)
      continue;

    if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx]
        && ComparisonTraitsT::less(_simplexTableau._reducedCosts[columnIdx], 0.0)) {
      enteringColumnIdx = columnIdx;
      spdlog::debug("ENTERING COLUMN IDX {} REDUCED COST {}", columnIdx, _simplexTableau._reducedCosts[columnIdx]);

      const std::vector<T> enteringColumn = computeEnteringColumn(*enteringColumnIdx);
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
RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::computeEnteringColumn(const int enteringColumnIdx) {
  std::vector<T> result(_simplexTableau._rowInfos.size());

  // TODO - maybe add kahan summation algo, maybe opt order
  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i)
  {
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
  const int columnBasicIdx = _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
  {
    if (j == columnBasicIdx)
    {
      result[j] = 1.0;
      continue;
    }

    if (j != columnBasicIdx && _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
    {
      result[j] = 0.0;
      continue;
    }

    for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k)
    {
      result[j] += _simplexTableau._basisMatrixInverse[rowIdx][k] *
                   _simplexTableau._constraintMatrix[k][j];
    }
  }

  return result;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::pivot(const int rowIdx, const int enteringColumnIdx,
                                                          const std::vector<T>& enteringColumn,
                                                          const std::vector<T>& pivotRow) {
  const PivotData<T> pivotData{rowIdx, enteringColumnIdx,
                               1.0 / enteringColumn[rowIdx]};
  spdlog::info("PIVOT COEFF {}, INV {}", enteringColumn[rowIdx], 1.0 / enteringColumn[rowIdx]);
  spdlog::info("RHS {}", _simplexTableau._rightHandSides[rowIdx]);
  updateReducedCosts(pivotData, pivotRow);
  updateInverseMatrixWithRHS(pivotData, enteringColumn);
  _simplexTableau.updateBasisData(pivotData);
}

template <typename T, typename ComparisonTraitsT>
std::optional<int> RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::chooseRowIdx(
    const std::vector<T>& enteringColumn) {
  std::optional<int> leavingRowIdx;

  // TODO - add choosing most negative

  for (int rowIdx = 0; rowIdx < _simplexTableau.getRowInfos().size(); ++rowIdx)
    if (ComparisonTraitsT::greater(enteringColumn[rowIdx], 0.0))
      if (!leavingRowIdx.has_value() ||
          ComparisonTraitsT::less(_simplexTableau._rightHandSides[rowIdx] *
                                      enteringColumn[*leavingRowIdx],
                                  _simplexTableau._rightHandSides[*leavingRowIdx] *
                                      enteringColumn[rowIdx]))
        leavingRowIdx = rowIdx;

  return leavingRowIdx;
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::updateReducedCosts(
    const PivotData<T> &pivotData, const std::vector<T> &pivotRow) {
  const auto& [leavingRowIdx, enteringColumnIdx, pivotingTermInverse] = pivotData;
  const auto commonCoeffReducedCost =
      _simplexTableau._reducedCosts[enteringColumnIdx] * pivotingTermInverse;
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    _simplexTableau._reducedCosts[j] -=
        commonCoeffReducedCost * pivotRow[j];
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::
    updateInverseMatrixWithRHS(const PivotData<T> &pivotData, const std::vector<T>& enteringColumn) {
  const auto& [leavingRowIdx, enteringColumnIdx, pivotingTermInverse] = pivotData;

  for (int i = 0; i < _simplexTableau._rowInfos.size(); ++i) {
    if (i == leavingRowIdx)
      continue;

    const auto commonCoeff = enteringColumn[i] * pivotingTermInverse;

    for (int j = 0; j < _simplexTableau._basisMatrixInverse.size(); ++j)
      _simplexTableau._basisMatrixInverse[i][j] -= commonCoeff * _simplexTableau._basisMatrixInverse[leavingRowIdx][j];

    _simplexTableau._rightHandSides[i] -= commonCoeff * _simplexTableau._rightHandSides[leavingRowIdx];
  }

  for (int j = 0; j < _simplexTableau._basisMatrixInverse.size(); ++j)
    _simplexTableau._basisMatrixInverse[leavingRowIdx][j] *= pivotingTermInverse;

  _simplexTableau._rightHandSides[leavingRowIdx] *= pivotingTermInverse;
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T,
                   ComparisonTraitsT>::removeArtificialVariablesFromProgram() {
  std::optional<int> firstArtificialIdx;
  for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    if (_simplexTableau._variableInfos[j]._isArtificial)
    {
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
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::removeArtificialVariablesFromBasis() {
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
  {
    const auto basicVarIdx = _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    spdlog::info("FOUND BASIC ARTIFICIAL COLUMN IDX {}", basicVarIdx);

    std::optional<int> nonZeroEntryColumnIndex;
    const std::vector<T> pivotRow = computePivotRow(rowIdx);

    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    {
      if (!_simplexTableau._variableInfos[j]._isArtificial
          || _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        continue;

      if (!ComparisonTraitsT::equal(pivotRow[j], 0.0))
      {
        nonZeroEntryColumnIndex = j;
        break;
      }
    }

    if (!nonZeroEntryColumnIndex.has_value())
    {
      spdlog::warn("Redundant constraints in lp formulation!");
//      removeRow(rowIdx);
    }
    else
    {
      const std::vector<T> enteringColumn = computeEnteringColumn(*nonZeroEntryColumnIndex);
      pivot(rowIdx, *nonZeroEntryColumnIndex, enteringColumn, pivotRow);
    }
  }
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::removeRow(const int rowIdx) {
  // rethink it
//  auto& [rowToBasisColumnIdxMap, isBasicColumnIndexBitset] = _simplexTableau._simplexBasisData;
//  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[rowIdx]] = false;
//  for (int i = rowIdx; i < _simplexTableau._rowInfos.size() - 1; ++i)
//    rowToBasisColumnIdxMap[i] = rowToBasisColumnIdxMap[i + 1];
//  rowToBasisColumnIdxMap.pop_back();
//
//  _simplexTableau._constraintMatrix.erase(_simplexTableau._constraintMatrix.begin() + rowIdx);
//  _simplexTableau._rowInfos.erase(_simplexTableau._rowInfos.begin() + rowIdx);
//  _simplexTableau._rightHandSides.erase(_simplexTableau._rightHandSides.begin() + rowIdx);
//  for (int i = 0; i < _simplexTableau._basisMatrixInverse.size(); ++i)
//    _simplexTableau._basisMatrixInverse[i].erase(_simplexTableau._basisMatrixInverse[i].begin() + rowIdx);
//  _simplexTableau._basisMatrixInverse.erase(_simplexTableau._basisMatrixInverse.begin() + rowIdx);
}
template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::setInitialObjective() {
  _simplexTableau._objectiveRow = _simplexTableau._initialProgram.getObjective();
  calculateDual();
  _simplexTableau.calculateReducedCosts();
  _simplexTableau.calculateCurrentObjectiveValue();
}

template <typename T, typename ComparisonTraitsT>
void RevisedPrimalSimplexPFI<T, ComparisonTraitsT>::calculateDual() {
  _simplexTableau._y.resize(_simplexTableau._basisMatrixInverse.size());
  for (int colIdx = 0; colIdx < _simplexTableau._y.size(); ++colIdx) {
    T sum{};
    for (int k = 0; k < _simplexTableau._rowInfos.size(); ++k)
      sum += _simplexTableau._objectiveRow[
_simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[k]] * _simplexTableau._basisMatrixInverse[k][colIdx];

    _simplexTableau._y[colIdx] = sum;
  }
}

template class RevisedPrimalSimplexPFI<double>;