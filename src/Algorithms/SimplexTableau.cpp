#include "SimplexTableau.h"

#include "src/DataModel/LinearProgram.h"
#include "src/Util/LPPrinter.h"
#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT>
SimplexTableau<T, SimplexTraitsT>::SimplexTableau(
    const LinearProgram<T> &linearProgram, const SimplexType simplexType,
    const SimplexTableauType simplexTableauType)
    : _initialProgram(linearProgram),
      _variableInfos(linearProgram._variableInfos),
      _variableLowerBounds(linearProgram._variableLowerBounds),
      _variableUpperBounds(linearProgram._variableUpperBounds),
      _variableLabelSet(linearProgram._variableLabelSet),
      _rowInfos(linearProgram._rowInfos),
      _constraintMatrix(linearProgram._constraintMatrix),
      _rightHandSides(linearProgram._rightHandSides),
      _initialRightHandSides(linearProgram._rightHandSides),
      _simplexTableauType(simplexTableauType),
      _result(LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
  if constexpr (SimplexTraitsT::useSparseRepresentationValue) {
    if (_simplexTableauType !=
        SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE) {
      SPDLOG_ERROR("SPARSE REPRESENTATION AVAILABLE ONLY FOR PFI");
    }
  }
}

template <typename T, typename SimplexTraitsT>
std::vector<T> SimplexTableau<T, SimplexTraitsT>::originalObjective() const {
  std::vector<T> result = _initialProgram.getObjective();
  result.resize(_variableInfos.size());
  return result;
}

template <typename T, typename SimplexTraitsT>
std::string SimplexTableau<T, SimplexTraitsT>::toString() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printLineBreak();

  if (!_simplexBasisData._rowToBasisColumnIdxMap.empty()) {
    lpPrinter.printVariableInfos(std::cref(_simplexBasisData));
    lpPrinter.printLineBreak();
  }

  if (!_reducedCosts.empty())
    lpPrinter.printReducedCostWithObjectiveValue(_reducedCosts,
                                                 _objectiveValue);

  if (_simplexTableauType == SimplexTableauType::FULL) {
    lpPrinter.printMatrixWithRHS(_simplexBasisData._rowToBasisColumnIdxMap,
                                 _fullTableau, _rightHandSides);
  } else {
    lpPrinter.printMatrixWithRHS1(_simplexBasisData._rowToBasisColumnIdxMap,
                                  _constraintMatrix, _rightHandSides,
                                  _initialRightHandSides);
  }
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

template <typename T, typename SimplexTraitsT>
std::string SimplexTableau<T, SimplexTraitsT>::toStringSolution() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printSolution(_x);
  return lpPrinter.toString();
}

template <typename T, typename SimplexTraitsT>
std::string SimplexTableau<T, SimplexTraitsT>::toStringSolutionWithDual(
    const LinearProgram<T> &dualProgram) const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printDual(_y, dualProgram);
  lpPrinter.printLineBreak();
  lpPrinter.printSolution(_x);
  lpPrinter.printLineBreak();
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

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateCurrentObjectiveValue() {
  typename SimplexTraitsT::CurrentAdder adder;
  for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx) {
    adder.addValue(_x[varIdx] * _objectiveRow[varIdx]);
  }

  _objectiveValue = adder.currentSum();
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateSolution() {
  _x.resize(_variableInfos.size());
  for (int colIdx = 0; colIdx < _variableInfos.size(); ++colIdx)
    if (!_simplexBasisData._isBasicColumnIndexBitset[colIdx])
      _x[colIdx] = *curSatisfiedBound(colIdx);

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx)
    _x[basicColumnIdx(rowIdx)] = _rightHandSides[rowIdx];
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::initTableau() {
  switch (_simplexTableauType) {
  case SimplexTableauType::FULL:
    return initFullTableau();
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE:
    return;
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return initBasisMatrixInverse();
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::initBasisMatrixInverse() {
  _basisMatrixInverse.resize(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _basisMatrixInverse.size(); ++rowIdx) {
    _basisMatrixInverse[rowIdx].resize(_basisMatrixInverse.size());
    _basisMatrixInverse[rowIdx][rowIdx] = 1;
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::initFullTableau() {
  _fullTableau = _constraintMatrix;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateDual() {
  switch (_simplexTableauType) {
  case SimplexTableauType::FULL:
    return calculateDualExplicit();
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue)
      return calculateDualPFISparse();

    return calculateDualPFI();
  }
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return calculateDualExplicit();
  }

  return calculateDualExplicit();
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateDualExplicit() {
  _y.resize(_basisMatrixInverse.size());
  for (int colIdx = 0; colIdx < _y.size(); ++colIdx) {
    typename SimplexTraitsT::CurrentAdder adder;

    for (int k = 0; k < _rowInfos.size(); ++k)
      adder.addValue(_objectiveRow[basicColumnIdx(k)] *
                     _basisMatrixInverse[k][colIdx]);

    _y[colIdx] = adder.currentSum();
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateDualPFI() {
  _y.resize(_rowInfos.size());
  for (int k = 0; k < _rowInfos.size(); ++k)
    _y[k] = _objectiveRow[basicColumnIdx(k)];

  multiplyByBasisMatrixRightInverseUsingPFI(_y);
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateDualPFISparse() {
  _y.resize(_rowInfos.size());
  for (int k = 0; k < _rowInfos.size(); ++k)
    _y[k] = _objectiveRow[basicColumnIdx(k)];

  multiplyByBasisMatrixRightInverseUsingPFISparseNormal(_y);
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
void SimplexTableau<T, SimplexTraitsT>::calculateReducedCostsFullTableau() {
  std::vector<T> objectiveBasisIndices(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx)
    objectiveBasisIndices[rowIdx] = _objectiveRow[basicColumnIdx(rowIdx)];

  _reducedCosts.resize(_variableInfos.size());
  for (int columnIdx = 0; columnIdx < _variableInfos.size(); ++columnIdx) {
    T yAn = _objectiveRow[columnIdx];
    for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx)
      yAn -= objectiveBasisIndices[rowIdx] * _fullTableau[rowIdx][columnIdx];

    _reducedCosts[columnIdx] = yAn;
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHS() {
  switch (_simplexTableauType) {
  case SimplexTableauType::FULL:
    return calculateRHSFullTableau();
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue)
      return calculateRHSPFISparse();

    return calculateRHSPFI();
  }
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return calculateRHSBasisInverse();
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHSFullTableau() {
  calculateRHSWithoutInverse();
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHSBasisInverse() {
  calculateRHSWithoutInverse();
  std::vector<T> tempRHS = _rightHandSides;
  for (int i = 0; i < _rowInfos.size(); ++i) {
    _rightHandSides[i] =
        SimplexTraitsT::dotProduct(_basisMatrixInverse[i], tempRHS);
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHSPFI() {
  calculateRHSWithoutInverse();
  multiplyByBasisMatrixLeftInverseUsingPFI(_rightHandSides);
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHSPFISparse() {
  calculateRHSWithoutInverse();
  multiplyByBasisMatrixLeftInverseUsingPFISparseNormal(_rightHandSides);
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHSWithoutInverse() {
  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    typename SimplexTraitsT::CurrentAdder adder;
    adder.addValue(_initialRightHandSides[rowIdx]);
    for (int j = 0; j < _variableInfos.size(); ++j)
      if (!_simplexBasisData._isBasicColumnIndexBitset[j])
        adder.addValue(-(*curSatisfiedBound(j) * _constraintMatrix[rowIdx][j]));

    _rightHandSides[rowIdx] = adder.currentSum();
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateBasisData(
    const PivotData<T> &pivotData) {
  const auto &[leavingRowIdx, enteringColumnIdx, _] = pivotData;
  auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset, _1, _2, _3] =
      _simplexBasisData;
  SPDLOG_DEBUG("LEAVING VARIABLE ROW IDX {} COLUMN IDX", leavingRowIdx,
               rowToBasisColumnIdxMap[leavingRowIdx]);

  isBasicColumnIndexBitset[enteringColumnIdx] = true;
  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[leavingRowIdx]] = false;
  rowToBasisColumnIdxMap[leavingRowIdx] = enteringColumnIdx;
}

template <typename T, typename SimplexTraitsT>
bool SimplexTableau<T, SimplexTraitsT>::isColumnEligibleToEnterBasis(
    const int colIdx) {
  return !_simplexBasisData._isBasicColumnIndexBitset[colIdx] &&
         (!_simplexBasisData._isColumnEligibleBitset.has_value() ||
          _simplexBasisData._isColumnEligibleBitset.value()[colIdx]);
}

template <typename T, typename SimplexTraitsT>
auto SimplexTableau<T, SimplexTraitsT>::computeTableauColumnGeneric(
    const int colIdx) -> VectorT {
  if constexpr (SimplexTraitsT::useSparseRepresentationValue)
    return computeTableauColumnPFISparse(colIdx);
  else {
    switch (_simplexTableauType) {
    case SimplexTableauType::FULL:
      return retrieveTableauColumn(colIdx);
    case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {

      return computeTableauColumnPFI(colIdx);
    }
    case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
      return computeTableauColumnExplicit(colIdx);
    }
    return computeTableauColumnExplicit(colIdx);
  }
}

template <typename T, typename SimplexTraitsT>
std::vector<T>
SimplexTableau<T, SimplexTraitsT>::retrieveTableauColumn(const int colIdx) {
  std::vector<T> result(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx)
    result[rowIdx] = _fullTableau[rowIdx][colIdx];

  return result;
}

template <typename T, typename SimplexTraitsT>
std::vector<T> SimplexTableau<T, SimplexTraitsT>::computeTableauColumnExplicit(
    const int colIdx) {
  std::vector<T> result(_rowInfos.size());

  for (int i = 0; i < _rowInfos.size(); ++i)
    result[i] = SimplexTraitsT::dotProduct(
        _basisMatrixInverse[i], _constraintMatrixSparseForm._columns[colIdx]);

  return result;
}
template <typename T, typename SimplexTraitsT>
std::vector<T>
SimplexTableau<T, SimplexTraitsT>::computeTableauColumnPFI(const int colIdx) {
  std::vector<T> tableauColumn = _constraintMatrixNormalForm._columns[colIdx];
  multiplyByBasisMatrixLeftInverseUsingPFI(tableauColumn);
  return tableauColumn;
}
template <typename T, typename SimplexTraitsT>
SparseVector<T>
SimplexTableau<T, SimplexTraitsT>::computeTableauColumnPFISparse(
    const int colIdx) {
  SparseVector<T> tableauColumn = _constraintMatrixSparseForm._columns[colIdx];
  multiplyByBasisMatrixLeftInverseUsingPFISparse(tableauColumn);
  return tableauColumn;
}

template <typename T, typename SimplexTraitsT>
auto SimplexTableau<T, SimplexTraitsT>::computeTableauRowGeneric(
    const int rowIdx) -> std::shared_ptr<VectorT> {
  if constexpr (SimplexTraitsT::useSparseRepresentationValue)
    return computeTableauRowPFISparse(rowIdx);
  else {
    switch (_simplexTableauType) {
    case SimplexTableauType::FULL:
      return retrieveTableauRow(rowIdx);
    case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {

      return computeTableauRowPFI(rowIdx);
    }
    case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
      return computeTableauRowExplicit(rowIdx);
    }

    return computeTableauRowExplicit(rowIdx);
  }
}

template <typename T, typename SimplexTraitsT>
std::shared_ptr<std::vector<T>>
SimplexTableau<T, SimplexTraitsT>::retrieveTableauRow(const int rowIdx) {
  return std::shared_ptr<std::vector<T>>(std::shared_ptr<std::vector<T>>(),
                                         &_fullTableau[rowIdx]);
}

template <typename T, typename SimplexTraitsT>
std::shared_ptr<std::vector<T>>
SimplexTableau<T, SimplexTraitsT>::computeTableauRowExplicit(const int rowIdx) {
  std::shared_ptr<std::vector<T>> result =
      std::make_shared<std::vector<T>>(_variableInfos.size());
  const int columnBasicIdx = basicColumnIdx(rowIdx);
  for (int j = 0; j < _variableInfos.size(); ++j) {
    if (j == columnBasicIdx)
      (*result)[j] = 1.0;
    else if (j != columnBasicIdx &&
             _simplexBasisData._isBasicColumnIndexBitset[j])
      (*result)[j] = 0.0;
    else
      (*result)[j] = SimplexTraitsT::dotProduct(
          _basisMatrixInverse[rowIdx], _constraintMatrixSparseForm._columns[j]);
  }

  return result;
}
template <typename T, typename SimplexTraitsT>
std::shared_ptr<std::vector<T>>
SimplexTableau<T, SimplexTraitsT>::computeTableauRowPFI(const int rowIdx) {
  std::vector<T> eRowIdx(_rowInfos.size());
  eRowIdx[rowIdx] = 1.0;
  multiplyByBasisMatrixRightInverseUsingPFI(eRowIdx);

  const int columnBasicIdx = basicColumnIdx(rowIdx);
  std::shared_ptr<std::vector<T>> result =
      std::make_shared<std::vector<T>>(_variableInfos.size());
  for (int j = 0; j < _variableInfos.size(); ++j) {
    if (j == columnBasicIdx)
      (*result)[j] = 1.0;
    else if (j != columnBasicIdx &&
             _simplexBasisData._isBasicColumnIndexBitset[j])
      (*result)[j] = 0.0;
    else
      (*result)[j] = SimplexTraitsT::dotProduct(
          eRowIdx, _constraintMatrixNormalForm._columns[j]);
  }

  return result;
}
template <typename T, typename SimplexTraitsT>
std::shared_ptr<SparseVector<T>>
SimplexTableau<T, SimplexTraitsT>::computeTableauRowPFISparse(
    const int rowIdx) {
  SparseVector<T> eRowIdx;
  eRowIdx._indexedValues.push_back({T{1.0}, rowIdx});
  eRowIdx._normalVec.resize(_rowInfos.size());
  eRowIdx._normalVec[rowIdx] = 1.0;
  multiplyByBasisMatrixRightInverseUsingPFISparse(eRowIdx);

  const int columnBasicIdx = basicColumnIdx(rowIdx);
  std::shared_ptr<SparseVector<T>> result = std::make_shared<SparseVector<T>>();
  result->_normalVec.resize(_variableInfos.size());
  for (int j = 0; j < _variableInfos.size(); ++j) {
    if (j == columnBasicIdx) {
      result->_normalVec[j] = 1.0;
      result->_indexedValues.push_back({1.0, j});
    } else if (!_simplexBasisData._isBasicColumnIndexBitset[j]) {
      const auto product =
          SimplexTraitsT::template dotProduct<PositiveNegativeAdderSafe>(
              eRowIdx, _constraintMatrixSparseForm._columns[j]);
      if (!NumericalTraitsT::isZero(product)) {
        result->_normalVec[j] = product;
        result->_indexedValues.push_back({result->_normalVec[j], j});
      }
    }
  }

  return result;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateReducedCostsGeneric(
    const PivotData<T> &pivotData, const VectorT &pivotRow) {
  const auto &[_, enteringColumnIdx, pivotingTermInverse] = pivotData;
  const auto commonCoeffReducedCost =
      _reducedCosts[enteringColumnIdx] * pivotingTermInverse;
  for (int j = 0; j < _variableInfos.size(); ++j)
    _reducedCosts[j] = NumericalTraitsT::add(
        _reducedCosts[j], -(commonCoeffReducedCost * pivotRow[j]));
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateDualGeneric(
    const PivotData<T> &pivotData, const VectorT &pivotRow) {
  if (_y.empty()) {
    return;
  }
  //  SPDLOG_INFO("Y SIZE {}, ROWS COUNT {}", _y.size(), _rowInfos.size());

  const auto &[_, enteringColumnIdx, pivotingTermInverse] = pivotData;
  const auto commonCoeffReducedCost =
      _reducedCosts[enteringColumnIdx] * pivotingTermInverse;
  for (int j = 0; j < _y.size(); ++j)
    _y[j] =
        NumericalTraitsT::add(_y[j], -(commonCoeffReducedCost * pivotRow[j]));
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateInverseMatrixWithRHSGeneric(
    const PivotData<T> &pivotData, const VectorT &enteringColumn) {
  const auto &[leavingRowIdx, _, pivotingTermInverse] = pivotData;

  ElementaryMatrixT etm{enteringColumn, pivotingTermInverse, leavingRowIdx};
  if constexpr (SimplexTraitsT::useSparseRepresentationValue) {
    SimplexTraitsT::template multiplyByETM<
        typename NumericalTraitsT::SafeNumericalAddOp>(etm, _rightHandSides);
    _sparsePfiEtms.push_back(std::move(etm));
  } else {
    switch (_simplexTableauType) {
    case SimplexTableauType::FULL: {
      SimplexTraitsT::multiplyByETM(etm, _fullTableau);
      if (!_basisMatrixInverse.empty()) {
        SimplexTraitsT::multiplyByETM(etm, _basisMatrixInverse);
      }
      SimplexTraitsT::multiplyByETM(etm, _rightHandSides);
      break;
    }
    case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
      SimplexTraitsT::multiplyByETM(etm, _rightHandSides);
      _pfiEtms.push_back(std::move(etm));
      break;
    }
    case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE: {
      SimplexTraitsT::multiplyByETM(etm, _basisMatrixInverse);
      SimplexTraitsT::multiplyByETM(etm, _rightHandSides);
      break;
    }
    }
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::pivotGeneric(
    const int rowIdx, const int enteringColumnIdx,
    const VectorT &enteringColumn, const VectorT &pivotRow) {
  const PivotData<T> pivotData{rowIdx, enteringColumnIdx,
                               1.0 / enteringColumn[rowIdx]};
  SPDLOG_DEBUG("PIVOT VALUE {},  INV {}", enteringColumn[rowIdx],
               1.0 / enteringColumn[rowIdx]);

  updateReducedCostsGeneric(pivotData, pivotRow);
  updateDualGeneric(pivotData, pivotRow);
  updateInverseMatrixWithRHSGeneric(pivotData, enteringColumn);
  updateBasisData(pivotData);
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::pivotImplicitBoundsGeneric(
    const int pivotRowIdx, const int enteringColumnIdx,
    const VectorT &enteringColumn, const VectorT &pivotRow,
    const bool leavingVarBecomesLowerBound) {
  const auto leavingBasicColumnIdx = basicColumnIdx(pivotRowIdx);
  SPDLOG_DEBUG(
      "PIVOT - ENTERING COLUMN IDX {}, ROW IDX {} LEAVING COLUMN IDX {}",
      enteringColumnIdx, pivotRowIdx, leavingBasicColumnIdx);

  if (_variableInfos[leavingBasicColumnIdx]._isFree) {
    SPDLOG_ERROR("FREE VAR {} SHOULDN'T LEAVE BASIS", leavingBasicColumnIdx);
  }

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    _rightHandSides[rowIdx] = NumericalTraitsT::add(
        _rightHandSides[rowIdx],
        *curSatisfiedBound(enteringColumnIdx) * enteringColumn[rowIdx]);
  }
  _rightHandSides[pivotRowIdx] = NumericalTraitsT::add(
      _rightHandSides[pivotRowIdx],
      -(leavingVarBecomesLowerBound
            ? *_variableLowerBounds[leavingBasicColumnIdx]
            : *_variableUpperBounds[leavingBasicColumnIdx]));

  pivotGeneric(pivotRowIdx, enteringColumnIdx, enteringColumn, pivotRow);

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
std::optional<T>
SimplexTableau<T, SimplexTraitsT>::curSatisfiedBound(const int varIdx) {
  if (_simplexBasisData._isBasicColumnIndexBitset[varIdx]) {
    SPDLOG_ERROR("BASIC VAR IDX {} DOESN'T HAVE SATISFIED BOUND", varIdx);
  }

  if (_simplexBasisData._isColumnAtLowerBoundBitset[varIdx])
    return *_variableLowerBounds[varIdx];

  if (_simplexBasisData._isColumnAtUpperBoundBitset[varIdx])
    return *_variableUpperBounds[varIdx];

  if (_variableInfos[varIdx]._isFree) {
    return 0.0;
  }

  return std::nullopt;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::initMatrixRepresentations() {
  _constraintMatrixNormalForm.clear();
  _constraintMatrixSparseForm.clear();

  _constraintMatrixSparseForm._rows.resize(_rowInfos.size());
  _constraintMatrixSparseForm._columns.resize(_variableInfos.size());

  _constraintMatrixNormalForm._rows.resize(_rowInfos.size());
  for (int i = 0; i < _rowInfos.size(); ++i) {
    _constraintMatrixNormalForm._rows[i].resize(_variableInfos.size());
    _constraintMatrixSparseForm._rows[i]._normalVec.resize(
        _variableInfos.size());
  }

  _constraintMatrixNormalForm._columns.resize(_variableInfos.size());
  for (int j = 0; j < _variableInfos.size(); ++j) {
    _constraintMatrixNormalForm._columns[j].resize(_rowInfos.size());
    _constraintMatrixSparseForm._columns[j]._normalVec.resize(_rowInfos.size());
  }

  for (int i = 0; i < _rowInfos.size(); ++i)
    for (int j = 0; j < _variableInfos.size(); ++j) {
      const auto currentValue = _constraintMatrix[i][j];
      _constraintMatrixNormalForm._rows[i][j] = currentValue;
      _constraintMatrixNormalForm._columns[j][i] = currentValue;

      if (currentValue != T{}) {
        _constraintMatrixSparseForm._rows[i]._normalVec[j] = currentValue;
        _constraintMatrixSparseForm._rows[i]._indexedValues.push_back(
            {currentValue, j});
        _constraintMatrixSparseForm._columns[j]._normalVec[i] = currentValue;
        _constraintMatrixSparseForm._columns[j]._indexedValues.push_back(
            {currentValue, i});
      }
    }
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::setObjective(
    const std::vector<T> &newObjective) {
  _objectiveRow = newObjective;
  _objectiveRow.resize(_variableInfos.size());
  updateTableauForNewObjective();
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateTableauForNewObjective() {
  switch (_simplexTableauType) {
  case SimplexTableauType::FULL: {
    calculateReducedCostsFullTableau();
    break;
  }
  default: {
    calculateDual();
    calculateReducedCostsBasedOnDual();
    break;
  }
  }
  calculateSolution();
  calculateCurrentObjectiveValue();
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::
    multiplyByBasisMatrixLeftInverseUsingPFI(std::vector<T> &vec) {
  for (const auto &pfiEtm : _pfiEtms)
    SimplexTraitsT::multiplyByETM(pfiEtm, vec);
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::
    multiplyByBasisMatrixLeftInverseUsingPFISparse(SparseVector<T> &vec) {
  for (const auto &pfiEtm : _sparsePfiEtms)
    SimplexTraitsT::template multiplyByETM<
        typename NumericalTraitsT::SafeNumericalAddOp>(pfiEtm, vec);
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::
    multiplyByBasisMatrixLeftInverseUsingPFISparseNormal(std::vector<T> &vec) {
  for (const auto &pfiEtm : _sparsePfiEtms)
    SimplexTraitsT::template multiplyByETM<
        typename NumericalTraitsT::SafeNumericalAddOp>(pfiEtm, vec);
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::
    multiplyByBasisMatrixRightInverseUsingPFI(std::vector<T> &vec) {
  for (int i = 0; i < _pfiEtms.size(); ++i)
    SimplexTraitsT::multiplyByETMFromRight(vec,
                                           _pfiEtms[_pfiEtms.size() - 1 - i]);
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::
    multiplyByBasisMatrixRightInverseUsingPFISparse(SparseVector<T> &vec) {
  for (int i = 0; i < _sparsePfiEtms.size(); ++i)
    SimplexTraitsT::template multiplyByETMFromRight<PositiveNegativeAdderSafe>(
        vec, _sparsePfiEtms[_sparsePfiEtms.size() - 1 - i]);
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::
    multiplyByBasisMatrixRightInverseUsingPFISparseNormal(std::vector<T> &vec) {
  for (int i = 0; i < _sparsePfiEtms.size(); ++i)
    SimplexTraitsT::template multiplyByETMFromRight<PositiveNegativeAdderSafe>(
        vec, _sparsePfiEtms[_sparsePfiEtms.size() - 1 - i]);
}

template <typename T, typename SimplexTraitsT>
boost::dynamic_bitset<>
SimplexTableau<T, SimplexTraitsT>::getIsVariableFreeBitset() const {
  boost::dynamic_bitset<> result(_variableInfos.size());
  for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx) {
    result[varIdx] = _variableInfos[varIdx]._isFree;
  }
  return result;
}

template class SimplexTableau<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableau<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class SimplexTableau<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableau<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
