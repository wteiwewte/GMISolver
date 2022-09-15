#include "SimplexTableau.h"

#include "src/DataModel/LinearProgram.h"
#include "src/Util/LPPrinter.h"
#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT>
SimplexTableau<T, SimplexTraitsT>::SimplexTableau(const LinearProgram<T> &linearProgram,
                                  const bool isPrimalSimplex)
    : _initialProgram(linearProgram),
      _variableInfos(linearProgram._variableInfos),
      _variableLowerBounds(linearProgram._variableLowerBounds),
      _variableUpperBounds(linearProgram._variableUpperBounds),
      _variableLabelSet(linearProgram._variableLabelSet),
      _rowInfos(linearProgram._rowInfos),
      _constraintMatrix(linearProgram._constraintMatrix),
      _rightHandSides(linearProgram._rightHandSides),
      _initialRightHandSides(linearProgram._rightHandSides),
      _initialObjectiveRow(linearProgram._objective),
      _result(LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
  SPDLOG_INFO("Converting LP to standard form");
  convertToStandardForm();
  SPDLOG_TRACE(toString());

  if (isPrimalSimplex) {
    SPDLOG_INFO("Making RHS non-negative");
    makeRightHandSidesNonNegative();
    SPDLOG_TRACE(toString());
  }

  SPDLOG_INFO("Adding artificial variables");
  addArtificialVariables();
  _objectiveRow = isPrimalSimplex ? initialPrimalSimplexObjective()
                                  : initialDualSimplexObjective();
  init(isPrimalSimplex);
  SPDLOG_TRACE("Simplex tableau with artificial variables");
  SPDLOG_TRACE(toString());
  SPDLOG_TRACE(toStringLpSolveFormat());
  initMatrixRepresentations();
}

template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::addArtificialVariables() {
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
        newArtificialLabelStr, VariableType::CONTINUOUS, false, true, false});
  }
}

template <typename T, typename SimplexTraitsT>
std::vector<T> SimplexTableau<T, SimplexTraitsT>::initialPrimalSimplexObjective() const {
  std::vector<T> result(_constraintMatrix.front().size());

  for (int variableIdx = 0; variableIdx < _variableInfos.size(); ++variableIdx)
    if (_variableInfos[variableIdx]._isArtificial)
      result[variableIdx] = 1;

  return result;
}

template <typename T, typename SimplexTraitsT>
std::vector<T> SimplexTableau<T, SimplexTraitsT>::initialDualSimplexObjective() const {
  std::vector<T> result = _initialProgram.getObjective();
  result.resize(_variableInfos.size());
  return result;
}

template <typename T, typename SimplexTraitsT>
std::optional<SimplexBasisData>
SimplexTableau<T, SimplexTraitsT>::createBasisFromArtificialVars() const {
  std::optional<int> firstArtificialIdx;
  for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx)
    if (_variableInfos[varIdx]._isArtificial) {
      firstArtificialIdx = varIdx;
      break;
    }

  if (!firstArtificialIdx.has_value()) {
    SPDLOG_WARN("No artificial variable found");
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

template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::init(const bool isPrimalSimplex) {
  if (auto simplexBasisData = createBasisFromArtificialVars();
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

template <typename T, typename SimplexTraitsT> std::string SimplexTableau<T, SimplexTraitsT>::toString() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printLineBreak();

  if (!_simplexBasisData._rowToBasisColumnIdxMap.empty()) {
    lpPrinter.printVariableInfos(std::cref(_simplexBasisData));
    lpPrinter.printLineBreak();
  }

  if (!_reducedCosts.empty())
    lpPrinter.printReducedCostWithObjectiveValue(_reducedCosts,
                                                 _objectiveValue);

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
template <typename T, typename SimplexTraitsT>
std::string SimplexTableau<T, SimplexTraitsT>::toStringObjectiveValue() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printCurrentObjectiveValue(_result, _objectiveValue);
  return lpPrinter.toString();
}
template <typename T, typename SimplexTraitsT> std::string SimplexTableau<T, SimplexTraitsT>::toStringSolution() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printSolution(_x);
  return lpPrinter.toString();
}
template <typename T, typename SimplexTraitsT>
std::string SimplexTableau<T, SimplexTraitsT>::toStringLpSolveFormat() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printInLpSolveFormat(_constraintMatrix, _objectiveRow,
                                 _rightHandSides, _variableLowerBounds,
                                 _variableUpperBounds);
  return lpPrinter.toString();
}

template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::calculateCurrentObjectiveValue() {
  _objectiveValue = T{};
  for (int i = 0; i < _rowInfos.size(); ++i)
    _objectiveValue +=
        _rightHandSides[i] *
        _objectiveRow[basicColumnIdx(i)];

  for (int j = 0; j < _variableInfos.size(); ++j)
    if (!_simplexBasisData._isBasicColumnIndexBitset[j])
      _objectiveValue += *curSatisfiedBound(j) *
                                         _objectiveRow[j];
}
template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::calculateSolution() {
  _x.resize(_variableInfos.size());
  for (int colIdx = 0; colIdx < _variableInfos.size(); ++colIdx)
    if (!_simplexBasisData._isBasicColumnIndexBitset[colIdx])
      _x[colIdx] = *curSatisfiedBound(colIdx);

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx)
    _x[basicColumnIdx(rowIdx)] =
        _rightHandSides[rowIdx];
}
template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::initBasisMatrixInverse() {
  _basisMatrixInverse.resize(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _basisMatrixInverse.size(); ++rowIdx) {
    _basisMatrixInverse[rowIdx].resize(_basisMatrixInverse.size());
    _basisMatrixInverse[rowIdx][rowIdx] = 1;
  }
}
template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::initDual() {
  _y.resize(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _y.size(); ++rowIdx) {
    const auto basicColumnIdx =
        _simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    _y[rowIdx] = _objectiveRow[basicColumnIdx];
  }
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateReducedCostsBasedOnDual() {
  _reducedCosts.resize(_variableInfos.size());
  for (int columnIdx = 0; columnIdx < _variableInfos.size(); ++columnIdx) {
    T yAn = _objectiveRow[columnIdx];
    for (int i = 0; i < _rowInfos.size(); ++i)
      yAn -= _y[i] * _constraintMatrix[i][columnIdx];

    _reducedCosts[columnIdx] = yAn;
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHS() {
  std::vector<T> tempRHS = _initialRightHandSides;
  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    for (int j = 0; j < _variableInfos.size(); ++j)
      if (!_simplexBasisData._isBasicColumnIndexBitset[j])
        tempRHS[rowIdx] -= *curSatisfiedBound(j) *
                           _constraintMatrix[rowIdx][j];
  }

  std::fill(_rightHandSides.begin(),
            _rightHandSides.end(), T{0.0});
  for (int i = 0; i < _rowInfos.size(); ++i)
    for (int j = 0; j < _rowInfos.size(); ++j)
      _rightHandSides[i] +=
          _basisMatrixInverse[i][j] * tempRHS[j];
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateBasisData(const PivotData<T> &pivotData) {
  const auto &[leavingRowIdx, enteringColumnIdx, _] = pivotData;
  auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset, _1, _2] =
      _simplexBasisData;
  SPDLOG_DEBUG("LEAVING VARIABLE ROW IDX {} COLUMN IDX", leavingRowIdx,
               rowToBasisColumnIdxMap[leavingRowIdx]);

  isBasicColumnIndexBitset[enteringColumnIdx] = true;
  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[leavingRowIdx]] = false;
  rowToBasisColumnIdxMap[leavingRowIdx] = enteringColumnIdx;
}
template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::initBoundsForDualSimplex() {
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

template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::convertToStandardForm() {
  int nonEqualityRows = 0;
  for (const auto &rowInfo : _rowInfos)
    if (!rowInfo.isEqualityRow())
      ++nonEqualityRows;

  const int variableCountAtTheStart = _variableInfos.size();
  const int newVariableCount = variableCountAtTheStart + nonEqualityRows;
  int newVariableIdx = variableCountAtTheStart;

  const auto newSlackLabel = [&]() {
    const std::string firstPattern = "S" + std::to_string(newVariableIdx);
    return (_variableLabelSet.find(firstPattern) == _variableLabelSet.end())
               ? firstPattern
               : firstPattern + Constants::SLACK_SUFFIX;
  };

  _initialObjectiveRow.resize(newVariableCount);
  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    _constraintMatrix[rowIdx].resize(newVariableCount);
    if (_rowInfos[rowIdx].isEqualityRow())
      continue;

    switch (_rowInfos[rowIdx]._type) {
    // TODO - check variable type integer or continuous
    // TODO - fix slack labeling
    case RowType::LESS_THAN_OR_EQUAL: {
      _constraintMatrix[rowIdx][newVariableIdx] = 1;
      break;
    }
    case RowType::GREATER_THAN_OR_EQUAL: {
      _constraintMatrix[rowIdx][newVariableIdx] = -1;
      break;
    }
    default:
      break;
    }
    const auto newSlackLabelStr = newSlackLabel();
    _variableInfos.push_back(
        VariableInfo{newSlackLabelStr, VariableType::CONTINUOUS, true, false, false});
    _variableLabelSet.insert(newSlackLabelStr);
    _variableLowerBounds.push_back(0.0);
    _variableUpperBounds.push_back(std::nullopt);
    ++newVariableIdx;
    _rowInfos[rowIdx]._type = RowType::EQUALITY;
  }
}

template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::makeRightHandSidesNonNegative() {
  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    T sum{0.0};
    for (int variableIdx = 0; variableIdx < _variableInfos.size();
         ++variableIdx)
      sum += _variableLowerBounds[variableIdx].value() *
             _constraintMatrix[rowIdx][variableIdx];

    const T diff = _rightHandSides[rowIdx] - sum;

    if (diff < 0.0) {
      for (int variableIdx = 0; variableIdx < _variableInfos.size();
           ++variableIdx)
        _constraintMatrix[rowIdx][variableIdx] =
            -_constraintMatrix[rowIdx][variableIdx];
      _rightHandSides[rowIdx] = -_rightHandSides[rowIdx];
      _initialRightHandSides[rowIdx] = -_initialRightHandSides[rowIdx];
    }
  }
}

template <typename T, typename SimplexTraitsT> void SimplexTableau<T, SimplexTraitsT>::addBoundsToMatrix() {
  for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx) {
    if (const auto lowerBound = _variableLowerBounds[varIdx];
        lowerBound.has_value()) {
      _constraintMatrix.emplace_back(_variableInfos.size())[varIdx] = 1;
      _rightHandSides.push_back(*lowerBound);
      _rowInfos.push_back(RowInfo{{}, RowType::GREATER_THAN_OR_EQUAL});
    }

    if (const auto upperBound = _variableUpperBounds[varIdx];
        upperBound.has_value()) {
      _constraintMatrix.emplace_back(_variableInfos.size())[varIdx] = 1;
      _rightHandSides.push_back(*upperBound);
      _rowInfos.push_back(RowInfo{{}, RowType::LESS_THAN_OR_EQUAL});
    }
  }
}

template <typename T, typename SimplexTraitsT>
bool SimplexTableau<T, SimplexTraitsT>::isColumnAllowedToEnterBasis(const int colIdx) {
  return  !_variableInfos[colIdx]._isArtificial &&
          !_variableInfos[colIdx]._isFixed &&
          !_simplexBasisData._isBasicColumnIndexBitset[colIdx];
}
template <typename T, typename SimplexTraitsT>
std::vector<T>
SimplexTableau<T, SimplexTraitsT>::computeTableauColumn(const int colIdx) {
  std::vector<T> result(_rowInfos.size());

  for (int i = 0; i < _rowInfos.size(); ++i)
    result[i] = SimplexTraitsT::dotProduct(_basisMatrixInverse[i], _constraintMatrixSparseForm._columns[colIdx]);

  return result;
}

template <typename T, typename SimplexTraitsT>
std::vector<T> SimplexTableau<T, SimplexTraitsT>::computeTableauRow(const int rowIdx) {
  std::vector<T> result(_variableInfos.size());

  // TODO - maybe add kahan summation algo, maybe opt order
  const int columnBasicIdx = basicColumnIdx(rowIdx);
  for (int j = 0; j < _variableInfos.size(); ++j) {
    if (j == columnBasicIdx)
      result[j] = 1.0;
    else if (j != columnBasicIdx &&
             _simplexBasisData._isBasicColumnIndexBitset[j])
      result[j] = 0.0;
    else
      result[j] = SimplexTraitsT::dotProduct(_basisMatrixInverse[rowIdx], _constraintMatrixSparseForm._columns[j]);
  }

  return result;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::pivot(const int rowIdx, const int enteringColumnIdx,
                              const std::vector<T> &enteringColumn,
                              const std::vector<T> &pivotRow) {
  const PivotData<T> pivotData{rowIdx, enteringColumnIdx,
                               1.0 / enteringColumn[rowIdx]};
  SPDLOG_DEBUG("PIVOT VALUE {},  INV {}", enteringColumn[rowIdx],
               1.0 / enteringColumn[rowIdx]);

  updateReducedCosts(pivotData, pivotRow);
  updateInverseMatrixWithRHS(pivotData, enteringColumn);
  updateBasisData(pivotData);
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::pivotImplicitBounds(
    const int pivotRowIdx, const int enteringColumnIdx,
    const std::vector<T> &enteringColumn, const std::vector<T> &pivotRow,
    const bool leavingVarBecomesLowerBound) {
  const auto leavingBasicColumnIdx = basicColumnIdx(pivotRowIdx);
  SPDLOG_DEBUG(
      "PIVOT - ENTERING COLUMN IDX {}, ROW IDX {} LEAVING COLUMN IDX {}",
      enteringColumnIdx, pivotRowIdx, leavingBasicColumnIdx);

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    _rightHandSides[rowIdx] +=
        *curSatisfiedBound(enteringColumnIdx) *
        enteringColumn[rowIdx];

  }
  _rightHandSides[pivotRowIdx] -=
      (leavingVarBecomesLowerBound
           ? *_variableLowerBounds[leavingBasicColumnIdx]
           : *_variableUpperBounds[leavingBasicColumnIdx]);

  pivot(pivotRowIdx, enteringColumnIdx, enteringColumn, pivotRow);

  auto &isColumnAtLowerBoundBitset =
      _simplexBasisData._isColumnAtLowerBoundBitset;
  auto &isColumnAtUpperBoundBitset =
      _simplexBasisData._isColumnAtUpperBoundBitset;

  (leavingVarBecomesLowerBound
       ? isColumnAtLowerBoundBitset[leavingBasicColumnIdx]
       : isColumnAtUpperBoundBitset[leavingBasicColumnIdx]) = true;

  if (isColumnAtLowerBoundBitset[enteringColumnIdx])
    isColumnAtLowerBoundBitset[enteringColumnIdx] = false;

  if (isColumnAtUpperBoundBitset[enteringColumnIdx])
    isColumnAtUpperBoundBitset[enteringColumnIdx] = false;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateReducedCosts(const PivotData<T> &pivotData,
                                           const std::vector<T> &pivotRow) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;
  const auto commonCoeffReducedCost =
      _reducedCosts[enteringColumnIdx] * pivotingTermInverse;
  for (int j = 0; j < _variableInfos.size(); ++j)
    _reducedCosts[j] -= commonCoeffReducedCost * pivotRow[j];
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateInverseMatrixWithRHS(
    const PivotData<T> &pivotData, const std::vector<T> &enteringColumn) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;

  const ElementaryMatrix<T> etm{enteringColumn,pivotingTermInverse, leavingRowIdx};

  SimplexTraitsT::multiplyByETM(etm, _basisMatrixInverse);
  SimplexTraitsT::multiplyByETM(etm, _rightHandSides);
}
template <typename T, typename SimplexTraitsT>
std::optional<T> SimplexTableau<T, SimplexTraitsT>::curSatisfiedBound(const int varIdx) {
  if (_simplexBasisData._isColumnAtLowerBoundBitset[varIdx])
    return *_variableLowerBounds[varIdx];

  if (_simplexBasisData._isColumnAtUpperBoundBitset[varIdx])
    return *_variableUpperBounds[varIdx];

  return std::nullopt;
}

template <typename T, typename SimplexTraitsT>
bool SimplexTableau<T, SimplexTraitsT>::reinversion() {
  // TODO - opt it by storing columns in vectors
  const auto basisSize = _rowInfos.size();
  std::vector<int> columnIndexesMapping(basisSize);
  std::vector<std::vector<T>> basisColumns(basisSize);
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    const auto basicColumnIdx =
        _simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
    columnIndexesMapping[rowIdx] = basicColumnIdx;
    basisColumns[rowIdx].resize(basisSize);
    for (int j = 0; j < basisSize; ++j)
      basisColumns[rowIdx][j] =
          _constraintMatrix[j][basicColumnIdx];
  }

  std::vector<std::vector<T>> newBasisMatrixInverse(basisSize);
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    newBasisMatrixInverse[rowIdx].resize(basisSize);
    newBasisMatrixInverse[rowIdx][rowIdx] = 1.0;
  }

  std::vector<bool> isUnusedColumn(basisSize, true);

  [[maybe_unused]] const auto findPivotColumnFirstEligible =
      [&](const int rowIdx) -> std::optional<int> {
    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx] &&
          !SimplexTraitsT::equal(basisColumns[colIdx][rowIdx], 0.0))
        return colIdx;

    return std::nullopt;
  };

  [[maybe_unused]] const auto findPivotColumnMaxAbsValue =
      [&](const int rowIdx) -> std::optional<int> {
    std::optional<T> maxAbsPivotValue;
    std::optional<int> maxAbsPivotColIdx;
    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx] &&
          SimplexTraitsT::isEligibleForPivot(basisColumns[colIdx][rowIdx])) {
        const auto currentAbsPivotValue =
            std::fabs(basisColumns[colIdx][rowIdx]);
        if (!maxAbsPivotColIdx.has_value() ||
            (*maxAbsPivotColIdx < currentAbsPivotValue)) {
          maxAbsPivotValue = currentAbsPivotValue;
          maxAbsPivotColIdx = colIdx;
        }
      }

    return maxAbsPivotColIdx;
  };

  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    const auto pivotColumnIdx = findPivotColumnMaxAbsValue(rowIdx);
    if (!pivotColumnIdx.has_value()) {
      SPDLOG_WARN("Basis matrix reinversion failed for row {}!", rowIdx);
      return false;
    }

    const T pivotingTermInverse{1.0 / basisColumns[*pivotColumnIdx][rowIdx]};
    for (int j = 0; j < basisSize; ++j) {
      if (j == rowIdx)
        continue;

      //      if (SimplexTraitsT::equal(basisColumns[*pivotColumnIdx][j],
      //      0.0))
      //        continue;

      const auto commonCoeff =
          pivotingTermInverse * basisColumns[*pivotColumnIdx][j];

      //      if (SimplexTraitsT::equal(commonCoeff, 0.0))
      //        continue;
      for (int k = 0; k < basisSize; ++k)
        newBasisMatrixInverse[j][k] -=
            commonCoeff * newBasisMatrixInverse[rowIdx][k];
    }

    for (int k = 0; k < basisSize; ++k)
      newBasisMatrixInverse[rowIdx][k] *= pivotingTermInverse;

    isUnusedColumn[*pivotColumnIdx] = false;
    _simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
        columnIndexesMapping[*pivotColumnIdx];

    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx]) {
        //        if (SimplexTraitsT::equal(basisColumns[colIdx][rowIdx],
        //        0.0))
        //          continue;

        const auto commonCoeff =
            pivotingTermInverse * basisColumns[colIdx][rowIdx];

        //        if (!SimplexTraitsT::equal(commonCoeff, 0.0))
        for (int k = 0; k < basisSize; ++k) {
          if (k == rowIdx)
            continue;

          basisColumns[colIdx][k] -=
              commonCoeff * basisColumns[*pivotColumnIdx][k];
        }

        basisColumns[colIdx][rowIdx] *= pivotingTermInverse;
      }
  }

  _basisMatrixInverse.swap(newBasisMatrixInverse);
  calculateRHS();
  SPDLOG_INFO("REINVERSION SUCCESS");
  return true;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::initMatrixRepresentations() {
  _constraintMatrixNormalForm.clear();
  _constraintMatrixSparseForm.clear();

  _constraintMatrixNormalForm._rows.resize(_rowInfos.size());
  for (int i = 0; i < _rowInfos.size(); ++i)
    _constraintMatrixNormalForm._rows[i].resize(_variableInfos.size());

  _constraintMatrixNormalForm._columns.resize(_variableInfos.size());
  for (int j = 0; j < _variableInfos.size(); ++j)
    _constraintMatrixNormalForm._columns[j].resize(_rowInfos.size());

  _constraintMatrixSparseForm._rows.resize(_rowInfos.size());
  _constraintMatrixSparseForm._columns.resize(_variableInfos.size());

  for (int i = 0; i < _rowInfos.size(); ++i)
    for (int j = 0; j < _variableInfos.size(); ++j)
    {
      const auto currentValue = _constraintMatrix[i][j];
      _constraintMatrixNormalForm._rows[i][j] = currentValue;
      _constraintMatrixNormalForm._columns[j][i] = currentValue;

      if (currentValue != T{})
      {
        _constraintMatrixSparseForm._rows[i]._elements.push_back({currentValue, j});
        _constraintMatrixSparseForm._columns[j]._elements.push_back({currentValue, i});
      }
    }

//  int rowNonzeros = 0;
//  for (int i = 0; i < _rowInfos.size(); ++i)
//    SPDLOG_INFO("{} ROW NONZEROS COUNT {}", i,
//                (rowNonzeros += _constraintMatrixSparseForm._rows[i]._elements.size(), _constraintMatrixSparseForm._rows[i]._elements.size()));
//
//  SPDLOG_INFO("TOTAL NONZEROS COUNT {}", rowNonzeros);
//
//  for (int j = 0; j < _variableInfos.size(); ++j)
//    SPDLOG_INFO("{} COLUMN NONZEROS COUNT {}", j, _constraintMatrixSparseForm._columns[j]._elements.size());
}

template class SimplexTableau<double>;
