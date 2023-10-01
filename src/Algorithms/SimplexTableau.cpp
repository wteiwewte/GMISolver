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
  SPDLOG_INFO("Converting LP to standard form");
  convertToStandardForm();
  SPDLOG_TRACE(toString());

  if (simplexType == SimplexType::PRIMAL) {
    SPDLOG_INFO("Making RHS non-negative");
    makeRightHandSidesNonNegative();
    SPDLOG_TRACE(toString());
  }

  SPDLOG_INFO("Adding artificial variables");
  addArtificialVariables();
  initMatrixRepresentations();
  init(simplexType);

  calculateRHS();
  calculateCurrentObjectiveValue();
  calculateSolution();
  SPDLOG_TRACE("Simplex tableau with artificial variables");
  SPDLOG_TRACE(toString());
  SPDLOG_TRACE(toStringLpSolveFormat());

  if constexpr (SimplexTraitsT::useSparseRepresentationValue) {
    if (_simplexTableauType !=
        SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE) {
      SPDLOG_ERROR("SPARSE REPRESENTATION AVAILABLE ONLY FOR PFI");
    }
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::addArtificialVariables() {
  const int variableCountAtTheStart = _variableInfos.size();
  const int newVariableCount = variableCountAtTheStart + _rowInfos.size() - 1;

  const auto newArtificialLabel = [&](const auto varIdx) {
    return "A" + std::to_string(varIdx + 1);
  };

  _constraintMatrix[0].resize(newVariableCount);
  for (int rowIdx = 1; rowIdx < _rowInfos.size(); ++rowIdx) {
    _constraintMatrix[rowIdx].resize(newVariableCount);

    const auto newVariableIdx = variableCountAtTheStart + rowIdx - 1;
    _constraintMatrix[rowIdx][newVariableIdx] = 1;
    const auto newArtificialLabelStr = newArtificialLabel(newVariableIdx);
    _variableLowerBounds.push_back(0.0);
    _variableUpperBounds.push_back(std::nullopt);
    _variableInfos.push_back(VariableInfo{
        newArtificialLabelStr, VariableType::CONTINUOUS, false, true, false});
  }
}

template <typename T, typename SimplexTraitsT>
std::vector<T> SimplexTableau<T, SimplexTraitsT>::artificialObjective() const {
  std::vector<T> result(_constraintMatrix.front().size());

  for (int variableIdx = 0; variableIdx < _variableInfos.size(); ++variableIdx)
    if (_variableInfos[variableIdx]._isArtificial)
      result[variableIdx] = 1;

  return result;
}

template <typename T, typename SimplexTraitsT>
std::vector<T> SimplexTableau<T, SimplexTraitsT>::originalObjective() const {
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
  result._rowToBasisColumnIdxMap.resize(_rowInfos.size());

  // objective row -> objective var
  result._rowToBasisColumnIdxMap[0] = 0;
  result._isBasicColumnIndexBitset[0] = true;
  result._isColumnAtLowerBoundBitset[0] = false;

  for (int rowIdx = 1; rowIdx < _rowInfos.size(); ++rowIdx) {
    const auto basicColumnIdx = *firstArtificialIdx + rowIdx - 1;
    result._rowToBasisColumnIdxMap[rowIdx] = basicColumnIdx;
    result._isBasicColumnIndexBitset[basicColumnIdx] = true;
    result._isColumnAtLowerBoundBitset[basicColumnIdx] = false;
  }

  return result;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::init(const SimplexType simplexType) {
  if (auto simplexBasisData = createBasisFromArtificialVars();
      simplexBasisData.has_value())
    _simplexBasisData = std::move(*simplexBasisData);

  initBasisMatrixInverse();
  setObjective((simplexType == SimplexType::PRIMAL) ? artificialObjective()
                                                    : originalObjective());

  if (simplexType == SimplexType::DUAL) {
    initBoundsForDualSimplex();
  }
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

  lpPrinter.printMatrixWithRHS1(_simplexBasisData._rowToBasisColumnIdxMap,
                                _constraintMatrix, _rightHandSides,
                                _initialRightHandSides);
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
  for (int i = 0; i < _rowInfos.size(); ++i)
    adder.addValue(_rightHandSides[i] * _objectiveRow[basicColumnIdx(i)]);

  for (int j = 0; j < _variableInfos.size(); ++j)
    if (!_simplexBasisData._isBasicColumnIndexBitset[j])
      adder.addValue(*curSatisfiedBound(j) * _objectiveRow[j]);

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
void SimplexTableau<T, SimplexTraitsT>::initBasisMatrixInverse() {
  _basisMatrixInverse.resize(_rowInfos.size());
  for (int rowIdx = 0; rowIdx < _basisMatrixInverse.size(); ++rowIdx) {
    _basisMatrixInverse[rowIdx].resize(_basisMatrixInverse.size());
    _basisMatrixInverse[rowIdx][rowIdx] = 1;
  }
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

  //  _reducedCosts.resize(_variableInfos.size()); // FIXME why this doesnt work
  //  for (int columnIdx = 0; columnIdx < _variableInfos.size(); ++columnIdx)
  //    _reducedCosts[columnIdx] = _objectiveRow[columnIdx]
  //                                                   -SimplexTraitsT::dotProduct(_y,
  //                                                   _constraintMatrixSparseForm._columns[columnIdx]);
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHS() {
  switch (_simplexTableauType) {
  case SimplexTableauType::FULL:
    return calculateRHSExplicit();
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue)
      return calculateRHSPFISparse();

    return calculateRHSPFI();
  }
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return calculateRHSExplicit();
  }

  return calculateRHSExplicit();
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::calculateRHSExplicit() {
  calculateRHSWithoutInverse();
  std::vector<T> tempRHS = _rightHandSides;
  for (int i = 0; i < _rowInfos.size(); ++i)
    _rightHandSides[i] =
        SimplexTraitsT::dotProduct(_basisMatrixInverse[i], tempRHS);
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
  auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset, _1, _2] =
      _simplexBasisData;
  SPDLOG_DEBUG("LEAVING VARIABLE ROW IDX {} COLUMN IDX", leavingRowIdx,
               rowToBasisColumnIdxMap[leavingRowIdx]);

  isBasicColumnIndexBitset[enteringColumnIdx] = true;
  isBasicColumnIndexBitset[rowToBasisColumnIdxMap[leavingRowIdx]] = false;
  rowToBasisColumnIdxMap[leavingRowIdx] = enteringColumnIdx;
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::initBoundsForDualSimplex() {
  for (int varIdx = 1; varIdx < _variableInfos.size(); ++varIdx) {
    if (_reducedCosts[varIdx] < 0.0) {
      _simplexBasisData._isColumnAtLowerBoundBitset[varIdx] = false;
      _simplexBasisData._isColumnAtUpperBoundBitset[varIdx] = true;
    }

    if (_variableInfos[varIdx]._isArtificial) {
      _variableUpperBounds[varIdx] = 0.0;
      _variableInfos[varIdx]._isFixed = true;
    }
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::convertToStandardForm() {
  int nonEqualityRows = 0;
  for (const auto &rowInfo : _rowInfos)
    if (!rowInfo.isEqualityRow())
      ++nonEqualityRows;

  const int variableCountAtTheStart = _variableInfos.size();
  const int newVariableCount = variableCountAtTheStart + nonEqualityRows;
  int newVariableIdx = variableCountAtTheStart;

  const auto newSlackLabel = [&]() {
    const std::string firstPattern = "S" + std::to_string(newVariableIdx + 1);
    return (_variableLabelSet.find(firstPattern) == _variableLabelSet.end())
               ? firstPattern
               : firstPattern + Constants::SLACK_SUFFIX;
  };

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
    _variableInfos.push_back(VariableInfo{
        newSlackLabelStr, VariableType::CONTINUOUS, true, false, false});
    _variableLabelSet.insert(newSlackLabelStr);
    _variableLowerBounds.push_back(0.0);
    _variableUpperBounds.push_back(std::nullopt);
    ++newVariableIdx;
    _rowInfos[rowIdx]._type = RowType::EQUALITY;
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::makeRightHandSidesNonNegative() {
  for (int rowIdx = 1; rowIdx < _rowInfos.size(); ++rowIdx) {
    typename SimplexTraitsT::CurrentAdder adder;
    for (int variableIdx = 1; variableIdx < _variableInfos.size();
         ++variableIdx)
      adder.addValue(_variableLowerBounds[variableIdx].value() *
                     _constraintMatrix[rowIdx][variableIdx]);

    const T diff =
        NumericalTraitsT::add(_rightHandSides[rowIdx], -adder.currentSum());

    if (diff < 0.0) {
      for (int variableIdx = 1; variableIdx < _variableInfos.size();
           ++variableIdx)
        _constraintMatrix[rowIdx][variableIdx] =
            -_constraintMatrix[rowIdx][variableIdx];
      _rightHandSides[rowIdx] = -_rightHandSides[rowIdx];
      _initialRightHandSides[rowIdx] = -_initialRightHandSides[rowIdx];
    }
  }
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::addBoundsToMatrix() {
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
bool SimplexTableau<T, SimplexTraitsT>::isColumnAllowedToEnterBasis(
    const int colIdx) {
  return !_variableInfos[colIdx]._isArtificial &&
         !_variableInfos[colIdx]._isFixed &&
         !_simplexBasisData._isBasicColumnIndexBitset[colIdx];
}

template <typename T, typename SimplexTraitsT>
auto SimplexTableau<T, SimplexTraitsT>::computeTableauColumnGeneric(
    const int colIdx) -> VectorT {
  if constexpr (SimplexTraitsT::useSparseRepresentationValue)
    return computeTableauColumnPFISparse(colIdx);
  else {
    switch (_simplexTableauType) {
    case SimplexTableauType::FULL:
      return computeTableauColumnExplicit(colIdx);
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
  //  SPDLOG_INFO("COLUMN NONZEROS BEFORE {} AFTER {}",
  //  _constraintMatrixSparseForm._columns[colIdx]._indexedValues.size(),
  //              tableauColumn._indexedValues.size());
  return tableauColumn;
}

template <typename T, typename SimplexTraitsT>
auto SimplexTableau<T, SimplexTraitsT>::computeTableauRowGeneric(
    const int rowIdx) -> VectorT {
  if constexpr (SimplexTraitsT::useSparseRepresentationValue)
    return computeTableauRowPFISparse(rowIdx);
  else {
    switch (_simplexTableauType) {
    case SimplexTableauType::FULL:
      return computeTableauRowExplicit(rowIdx);
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
std::vector<T>
SimplexTableau<T, SimplexTraitsT>::computeTableauRowExplicit(const int rowIdx) {
  std::vector<T> result(_variableInfos.size());
  const int columnBasicIdx = basicColumnIdx(rowIdx);
  for (int j = 0; j < _variableInfos.size(); ++j) {
    if (j == columnBasicIdx)
      result[j] = 1.0;
    else if (j != columnBasicIdx &&
             _simplexBasisData._isBasicColumnIndexBitset[j])
      result[j] = 0.0;
    else
      result[j] = SimplexTraitsT::dotProduct(
          _basisMatrixInverse[rowIdx], _constraintMatrixSparseForm._columns[j]);
  }

  return result;
}
template <typename T, typename SimplexTraitsT>
std::vector<T>
SimplexTableau<T, SimplexTraitsT>::computeTableauRowPFI(const int rowIdx) {
  std::vector<T> eRowIdx(_rowInfos.size());
  eRowIdx[rowIdx] = 1.0;
  multiplyByBasisMatrixRightInverseUsingPFI(eRowIdx);

  const int columnBasicIdx = basicColumnIdx(rowIdx);
  std::vector<T> result(_variableInfos.size());
  for (int j = 0; j < _variableInfos.size(); ++j) {
    if (j == columnBasicIdx)
      result[j] = 1.0;
    else if (j != columnBasicIdx &&
             _simplexBasisData._isBasicColumnIndexBitset[j])
      result[j] = 0.0;
    else
      result[j] = SimplexTraitsT::dotProduct(
          eRowIdx, _constraintMatrixNormalForm._columns[j]);
  }

  return result;
}
template <typename T, typename SimplexTraitsT>
SparseVector<T> SimplexTableau<T, SimplexTraitsT>::computeTableauRowPFISparse(
    const int rowIdx) {
  SparseVector<T> eRowIdx;
  eRowIdx._indexedValues.push_back({T{1.0}, rowIdx});
  eRowIdx._normalVec.resize(_rowInfos.size());
  eRowIdx._normalVec[rowIdx] = 1.0;
  multiplyByBasisMatrixRightInverseUsingPFISparse(eRowIdx);

  const int columnBasicIdx = basicColumnIdx(rowIdx);
  SparseVector<T> result;
  result._normalVec.resize(_variableInfos.size());
  for (int j = 0; j < _variableInfos.size(); ++j) {
    if (j == columnBasicIdx) {
      result._normalVec[j] = 1.0;
      result._indexedValues.push_back({1.0, j});
    } else if (!_simplexBasisData._isBasicColumnIndexBitset[j]) {
      //      const auto product = SimplexTraitsT::dotProductSparse(eRowIdx,
      //      _constraintMatrixSparseForm._columns[j]); if (std::fabs(product) <
      //      1e-10 && std::fabs(product) > 0.0)
      //      {
      //        SPDLOG_INFO("DOT PRODUCT SPARSE {} {} {} {} {}", product,
      //                    SimplexTraitsT::template
      //                    dotProductSparse<SimpleAdderSafe>(eRowIdx,
      //                    _constraintMatrixSparseForm._columns[j]),
      //                    SimplexTraitsT::template
      //                    dotProductSparse<PositiveNegativeAdderSafe>(eRowIdx,
      //                    _constraintMatrixSparseForm._columns[j]),
      //                    SimplexTraitsT::template
      //                    dotProductSparse<KahanAdderSafe>(eRowIdx,
      //                    _constraintMatrixSparseForm._columns[j]),
      //                    SimplexTraitsT::template
      //                    dotProductSparse<KahanAdderNormal>(eRowIdx,
      //                    _constraintMatrixSparseForm._columns[j])
      //                    );
      //
      //      }
      const auto product =
          SimplexTraitsT::template dotProduct<PositiveNegativeAdderSafe>(
              eRowIdx, _constraintMatrixSparseForm._columns[j]);
      if (!NumericalTraitsT::isZero(product)) {
        result._normalVec[j] = product;
        result._indexedValues.push_back({result._normalVec[j], j});
      }
    }
  }

  //  SPDLOG_INFO("ROW NONZEROS BEFORE {} AFTER {}",
  //  _constraintMatrixSparseForm._rows[rowIdx]._indexedValues.size(),
  //              result._indexedValues.size());

  return result;
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateReducedCostsGeneric(
    const PivotData<T> &pivotData, const VectorT &pivotRow) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;
  const auto commonCoeffReducedCost =
      _reducedCosts[enteringColumnIdx] * pivotingTermInverse;
  for (int j = 0; j < _variableInfos.size(); ++j)
    _reducedCosts[j] = NumericalTraitsT::add(
        _reducedCosts[j], -(commonCoeffReducedCost * pivotRow[j]));
}

template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::updateInverseMatrixWithRHSGeneric(
    const PivotData<T> &pivotData, const VectorT &enteringColumn) {
  const auto &[leavingRowIdx, enteringColumnIdx, pivotingTermInverse] =
      pivotData;

  ElementaryMatrixT etm{enteringColumn, pivotingTermInverse, leavingRowIdx};

  if constexpr (SimplexTraitsT::useSparseRepresentationValue) {
    _sparsePfiEtms.push_back(std::move(etm));
    SimplexTraitsT::template multiplyByETM<
        typename NumericalTraitsT::SafeNumericalAddOp>(_sparsePfiEtms.back(),
                                                       _rightHandSides);
  } else {
    if (_simplexTableauType !=
        SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE) {
      SimplexTraitsT::multiplyByETM(etm, _basisMatrixInverse);
      SimplexTraitsT::multiplyByETM(etm, _rightHandSides);
    } else {
      SimplexTraitsT::multiplyByETM(etm, _rightHandSides);
      _pfiEtms.push_back(std::move(etm));
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
  if (_simplexBasisData._isColumnAtLowerBoundBitset[varIdx])
    return *_variableLowerBounds[varIdx];

  if (_simplexBasisData._isColumnAtUpperBoundBitset[varIdx])
    return *_variableUpperBounds[varIdx];

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

  //  int rowNonzeros = 0;
  //  for (int i = 0; i < _rowInfos.size(); ++i)
  //    SPDLOG_INFO("{} ROW NONZEROS COUNT {}", i,
  //                (rowNonzeros +=
  //                _constraintMatrixSparseForm._rows[i]._indexedValues.size(),
  //                _constraintMatrixSparseForm._rows[i]._indexedValues.size()));
  //
  //  SPDLOG_INFO("TOTAL NONZEROS COUNT {}", rowNonzeros);
  //
  //  for (int j = 0; j < _variableInfos.size(); ++j)
  //    SPDLOG_INFO("{} COLUMN NONZEROS COUNT {}", j,
  //    _constraintMatrixSparseForm._columns[j]._indexedValues.size());
}
template <typename T, typename SimplexTraitsT>
void SimplexTableau<T, SimplexTraitsT>::setObjective(
    const std::vector<T> &newObjective) {
  _objectiveRow = newObjective;
  _objectiveRow.resize(_variableInfos.size());
  calculateDual();
  calculateReducedCostsBasedOnDual();
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

template class SimplexTableau<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableau<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class SimplexTableau<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class SimplexTableau<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
