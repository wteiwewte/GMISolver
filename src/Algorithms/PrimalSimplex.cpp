#include "src/Algorithms/PrimalSimplex.h"

#include "src/Algorithms/SimplexTableau.h"

template <typename T, typename ComparisonTraitsT>
PrimalSimplex<T, ComparisonTraitsT>::PrimalSimplex(
    SimplexTableau<T> &simplexTableau) : _simplexTableau(simplexTableau) {}

template <typename T, typename ComparisonTraitsT>
void PrimalSimplex<T, ComparisonTraitsT>::runPhaseOne() {
  run();
  removeArtificialVariablesFromBasis();
  removeArtificialVariablesFromProgram();

  if (ComparisonTraitsT::greater(_simplexTableau._objectiveValue, 0.0))
  {
    spdlog::info("Program with artificial variable has optimum greater than 0 - "
                 "initial program is infeasible");
    return;
  }
//
//  setInitialObjective();
//  run();
}

template <typename T, typename ComparisonTraitsT>
void PrimalSimplex<T, ComparisonTraitsT>::run() {
    int iterCount = 1;
    while (true) {
      spdlog::info("ITERATION {}", iterCount++);
      const bool iterResult = runOneIteration();
      if (iterResult)
        break;

      _simplexTableau.calculateCurrentObjectiveValue();
      _simplexTableau.calculateSolution();
//      spdlog::info("{}\n", _simplexTableau.toStringShort());
    }
}
template <typename T, typename ComparisonTraitsT>
bool PrimalSimplex<T, ComparisonTraitsT>::runOneIteration() {
  std::optional<int> enteringColumnIdx;
  for (int columnIdx = 0; columnIdx < _simplexTableau.getVariableInfos().size(); ++columnIdx) {
    if(_simplexTableau._variableInfos[columnIdx]._isArtificial)
      continue;

    if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[columnIdx]
            && ComparisonTraitsT::less(_simplexTableau._reducedCosts[columnIdx], 0.0)) {
      enteringColumnIdx = columnIdx;
      spdlog::debug("ENTERING COLUMN IDX {} REDUCED COST {}", columnIdx, _simplexTableau._reducedCosts[columnIdx]);


      const std::optional<int> rowIdx = chooseRowIdx(*enteringColumnIdx);
      if (!rowIdx.has_value()) {
        _simplexTableau._result = LPOptimizationResult::UNBOUNDED;
        return true;
      }

      pivot(*rowIdx, *enteringColumnIdx);
      return false;
    }
  }

  _simplexTableau._result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
  return true;
}

template <typename T, typename ComparisonTraitsT>
void PrimalSimplex<T, ComparisonTraitsT>::pivot(const int rowIdx, const int enteringColumnIdx) {
  const PivotData<T> pivotData{rowIdx, enteringColumnIdx,
                     1.0 / _simplexTableau._constraintMatrix[rowIdx][enteringColumnIdx]};
  spdlog::info("PIVOT COEFF {}, INV {}", _simplexTableau._constraintMatrix[rowIdx][enteringColumnIdx], 1.0 / _simplexTableau._constraintMatrix[rowIdx][enteringColumnIdx]);
  spdlog::info("RHS {}", _simplexTableau._rightHandSides[rowIdx]);
  _simplexTableau.updateReducedCosts(pivotData);
  _simplexTableau.updateConstraintMatrixWithRHS(pivotData);
  _simplexTableau.updateBasisData(pivotData);
}

template <typename T, typename ComparisonTraitsT>
std::optional<int> PrimalSimplex<T, ComparisonTraitsT>::chooseRowIdx(
    const int enteringColumnIdx) {
    std::optional<int> leavingRowIdx;

    for (int rowIdx = 0; rowIdx < _simplexTableau.getRowInfos().size(); ++rowIdx)
      if (ComparisonTraitsT::greater(_simplexTableau._constraintMatrix[rowIdx][enteringColumnIdx], 0.0))
        if (!leavingRowIdx.has_value() ||
            ComparisonTraitsT::less(_simplexTableau._rightHandSides[rowIdx] *
                    _simplexTableau._constraintMatrix[*leavingRowIdx][enteringColumnIdx],
                _simplexTableau._rightHandSides[*leavingRowIdx] *
                    _simplexTableau._constraintMatrix[rowIdx][enteringColumnIdx]))
          leavingRowIdx = rowIdx;

//    spdlog::debug("PIVOT ROW IDX {}"

    return leavingRowIdx;
}

template <typename T, typename ComparisonTraitsT>
void PrimalSimplex<T,
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
void PrimalSimplex<T, ComparisonTraitsT>::removeArtificialVariablesFromBasis() {
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx)
  {
    const auto basicVarIdx = _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    if (!_simplexTableau._variableInfos[basicVarIdx]._isArtificial)
      continue;

    std::optional<int> nonZeroEntryColumnIndex;
    for (int j = 0; j < _simplexTableau._variableInfos.size(); ++j)
    {
      if (!_simplexTableau._variableInfos[j]._isArtificial
          || _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[j])
        continue;

      if (!ComparisonTraitsT::equal(_simplexTableau._constraintMatrix[rowIdx][j], 0.0))
      {
        nonZeroEntryColumnIndex = j;
        break;
      }
    }

    if (!nonZeroEntryColumnIndex.has_value())
    {
      spdlog::warn("Redundant constraints in lp formulation!");
      removeRow(rowIdx);
    }

    pivot(rowIdx, *nonZeroEntryColumnIndex);
  }
}
template <typename T, typename ComparisonTraitsT>
void PrimalSimplex<T, ComparisonTraitsT>::removeRow(const int rowIdx) {
  auto& [rowToBasisColumnIdxMap, isBasicColumnIndexBitset] = _simplexTableau._simplexBasisData;
  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[rowIdx]] = false;
  for (int i = rowIdx; i < _simplexTableau._rowInfos.size() - 1; ++i)
    rowToBasisColumnIdxMap[i] = rowToBasisColumnIdxMap[i + 1];
  rowToBasisColumnIdxMap.pop_back();

  _simplexTableau._constraintMatrix.erase(_simplexTableau._constraintMatrix.begin() + rowIdx);
  _simplexTableau._rowInfos.erase(_simplexTableau._rowInfos.begin() + rowIdx);
  _simplexTableau._rightHandSides.erase(_simplexTableau._rightHandSides.begin() + rowIdx);
}
template <typename T, typename ComparisonTraitsT>
void PrimalSimplex<T, ComparisonTraitsT>::setInitialObjective() {
  _simplexTableau._objectiveRow = _simplexTableau._initialProgram.getObjective();
  _simplexTableau.initDual();
  _simplexTableau.calculateReducedCosts();
  _simplexTableau.calculateCurrentObjectiveValue();
}

template class PrimalSimplex<double>;