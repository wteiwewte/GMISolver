#include "SimplexTableau.h"

#include "src/DataModel/LinearProgram.h"
#include "src/Util/ComparisonTraits.h"
#include "src/Util/LPPrinter.h"

template <typename T>
SimplexTableau<T>::SimplexTableau(const LinearProgram<T> &linearProgram,
                                  const bool isPrimalSimplex)
    : _initialProgram(linearProgram),
      _variableInfos(linearProgram._variableInfos),
      _variableLowerBounds(linearProgram._variableLowerBounds),
      _variableUpperBounds(linearProgram._variableUpperBounds),
      _rowInfos(linearProgram._rowInfos),
      _constraintMatrix(linearProgram._constraintMatrix),
      _rightHandSides(linearProgram._rightHandSides),
      _initialRightHandSides(linearProgram._rightHandSides),
      _objectiveRow((addArtificialVariables(),
                     isPrimalSimplex ? initialPrimalSimplexObjective()
                                     : initialDualSimplexObjective())),
      _result(LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
  init(isPrimalSimplex);
}

template <typename T> void SimplexTableau<T>::addArtificialVariables() {
  const int variableCountAtTheStart = _variableInfos.size();
  const int newVariableCount = variableCountAtTheStart + _rowInfos.size();

  const auto newArtificialLabel = [&](const auto varIdx) {
    return "A" + std::to_string(varIdx);
  };

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    _constraintMatrix[rowIdx].resize(newVariableCount);

    const auto newVariableIdx = variableCountAtTheStart + rowIdx;
    _constraintMatrix[rowIdx][newVariableIdx] = 1;
    const auto newArtificialLabelStr = newArtificialLabel(newVariableIdx);
    _variableLowerBounds.push_back(0.0);
    _variableUpperBounds.push_back(std::nullopt);
    _variableInfos.push_back(VariableInfo{
        newArtificialLabelStr, VariableType::CONTINUOUS, false, true});
  }
}

template <typename T>
std::vector<T> SimplexTableau<T>::initialPrimalSimplexObjective() const {
  std::vector<T> result(_constraintMatrix.front().size());

  for (int variableIdx = 0; variableIdx < _variableInfos.size(); ++variableIdx)
    if (_variableInfos[variableIdx]._isArtificial)
      result[variableIdx] = 1;

  return result;
}

template <typename T>
std::vector<T> SimplexTableau<T>::initialDualSimplexObjective() const {
  std::vector<T> result = _initialProgram.getObjective();
  result.resize(_variableInfos.size());
  return result;
}

template <typename T>
std::optional<SimplexBasisData>
SimplexTableau<T>::createBasisFromArtificialVars(
    const bool isPrimalSimplex) const {
  std::optional<int> firstArtificialIdx;
  for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx)
    if (_variableInfos[varIdx]._isArtificial) {
      firstArtificialIdx = varIdx;
      break;
    }

  if (!firstArtificialIdx.has_value()) {
    spdlog::warn("No artificial variable found");
    return std::nullopt;
  }

  SimplexBasisData result;
  result._isBasicColumnIndexBitset.resize(_variableInfos.size());
  result._isColumnAtLowerBoundBitset.resize(_variableInfos.size(), true);
  result._isColumnAtUpperBoundBitset.resize(_variableInfos.size(), false);
  result._rowToBasisColumnIdxMap.resize(_variableInfos.size());

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    const auto basicColumnIdx = *firstArtificialIdx + rowIdx;
    result._rowToBasisColumnIdxMap[rowIdx] = basicColumnIdx;
    result._isBasicColumnIndexBitset[basicColumnIdx] = true;
    result._isColumnAtLowerBoundBitset[basicColumnIdx] = false;
  }

  return result;
}

template <typename T> void SimplexTableau<T>::init(const bool isPrimalSimplex) {
  if (auto simplexBasisData = createBasisFromArtificialVars(isPrimalSimplex);
      simplexBasisData.has_value())
    _simplexBasisData = std::move(*simplexBasisData);

  initDual();
  initBasisMatrixInverse();
  calculateReducedCostsBasedOnDual();
  if (!isPrimalSimplex)
    initBoundsForDualSimplex();
  calculateCurrentObjectiveValue();
  calculateSolution();
}

template <typename T> std::string SimplexTableau<T>::toString() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printLineBreak();
  lpPrinter.printVariableInfos(std::cref(_simplexBasisData));
  lpPrinter.printLineBreak();
  lpPrinter.printReducedCostWithObjectiveValue(_reducedCosts, _objectiveValue);
  lpPrinter.printMatrixWithRHS(_simplexBasisData._rowToBasisColumnIdxMap,
                               _constraintMatrix, _rightHandSides);
  lpPrinter.printVariableBounds(_variableLowerBounds, _variableUpperBounds);
  lpPrinter.printLineBreak();
  lpPrinter.printInverseBasis(_basisMatrixInverse);
  lpPrinter.printLineBreak();
  lpPrinter.printDual(_y);
  lpPrinter.printLineBreak();
  lpPrinter.printSolution(_x);
  lpPrinter.printLineBreak();

  return lpPrinter.toString();
}
template <typename T> std::string SimplexTableau<T>::toStringShort() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  //  lpPrinter.printSolution(_x);
  lpPrinter.printCurrentObjectiveValue(_result, _objectiveValue);
  return lpPrinter.toString();
}
template <typename T>
std::string SimplexTableau<T>::toStringShortWithSolution() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printSolution(_x);
  lpPrinter.printCurrentObjectiveValue(_result, _objectiveValue);
  return lpPrinter.toString();
}
template <typename T>
std::string SimplexTableau<T>::toStringLpSolveFormat() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printInLpSolveFormat(_constraintMatrix, _objectiveRow,
                                 _rightHandSides, _variableLowerBounds,
                                 _variableUpperBounds);
  return lpPrinter.toString();
}

template <typename T> void SimplexTableau<T>::calculateCurrentObjectiveValue() {
  _objectiveValue = T{};
  for (int i = 0; i < _rowInfos.size(); ++i)
    _objectiveValue +=
        _rightHandSides[i] *
        _objectiveRow[_simplexBasisData._rowToBasisColumnIdxMap[i]];
}
template <typename T> void SimplexTableau<T>::calculateSolution() {
  _x.resize(_variableInfos.size());
  std::fill(_x.begin(), _x.end(), T{});
  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx)
    _x[_simplexBasisData._rowToBasisColumnIdxMap[rowIdx]] =
        _rightHandSides[rowIdx];
}
template <typename T> void SimplexTableau<T>::initBasisMatrixInverse() {
  _basisMatrixInverse.resize(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _basisMatrixInverse.size(); ++rowIdx) {
    _basisMatrixInverse[rowIdx].resize(_basisMatrixInverse.size());
    _basisMatrixInverse[rowIdx][rowIdx] = 1;
  }
}
template <typename T> void SimplexTableau<T>::initDual() {
  _y.resize(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _y.size(); ++rowIdx) {
    const auto basicColumnIdx =
        _simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    _y[rowIdx] = _objectiveRow[basicColumnIdx];
  }
}
template <typename T>
void SimplexTableau<T>::calculateReducedCostsBasedOnDual() {
  _reducedCosts.resize(_variableInfos.size());
  for (int columnIdx = 0; columnIdx < _variableInfos.size(); ++columnIdx) {
    T yAn = _objectiveRow[columnIdx];
    for (int i = 0; i < _rowInfos.size(); ++i)
      yAn -= _y[i] * _constraintMatrix[i][columnIdx];

    _reducedCosts[columnIdx] = yAn;
  }
}

template <typename T>
void SimplexTableau<T>::updateBasisData(const PivotData<T> &pivotData) {
  const auto &[leavingRowIdx, enteringColumnIdx, _] = pivotData;
  auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset, _1, _2] =
      _simplexBasisData;
  spdlog::debug("LEAVING VARIABLE ROW IDX {} COLUMN IDX", leavingRowIdx,
                rowToBasisColumnIdxMap[leavingRowIdx]);

  isBasicColumnIndexBitset[enteringColumnIdx] = true;
  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[leavingRowIdx]] = false;
  rowToBasisColumnIdxMap[leavingRowIdx] = enteringColumnIdx;
}
template <typename T> void SimplexTableau<T>::initBoundsForDualSimplex() {
  for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx) {
    if (_reducedCosts[varIdx] < 0.0) {
      _simplexBasisData._isColumnAtLowerBoundBitset[varIdx] = false;
      _simplexBasisData._isColumnAtUpperBoundBitset[varIdx] = true;
    }

    if (_variableInfos[varIdx]._isArtificial) {
      _variableUpperBounds[varIdx] = 0.0;
    }
  }
}

template class SimplexTableau<double>;
