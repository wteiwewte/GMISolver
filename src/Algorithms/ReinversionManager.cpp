#include "src/Algorithms/ReinversionManager.h"

#include "src/Algorithms/SimplexTableau.h"

template <typename T, typename SimplexTraitsT>
ReinversionManager<T, SimplexTraitsT>::ReinversionManager(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    const int32_t reinversionFrequency)
    : _simplexTableau(simplexTableau),
      _reinversionFrequency(reinversionFrequency) {}

template <typename T, typename SimplexTraitsT>
bool ReinversionManager<T, SimplexTraitsT>::tryReinverse() {
  ++_iterCount;
  if (_reinversionFrequency && (_iterCount % _reinversionFrequency == 0)) {
    return reinverse();
  }
  return true;
}
template <typename T, typename SimplexTraitsT>
bool ReinversionManager<T, SimplexTraitsT>::reinverse() {
  const std::optional<bool> basisReinversionResult = reinverseBasis();
  if (!basisReinversionResult.has_value())
    return true;

  if (!basisReinversionResult.value())
    return basisReinversionResult.value();

  updateTableau();
  return true;
}

template <typename T, typename SimplexTraitsT>
std::optional<bool> ReinversionManager<T, SimplexTraitsT>::reinverseBasis() {
  switch (_simplexTableau._simplexTableauType) {
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue)
      return reinverseBasisPFISparse();

    return reinverseBasisPFI();
  }
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return reinverseBasisExplicit();
  default:
    return std::nullopt;
  }
}

template <typename T, typename SimplexTraitsT>
void ReinversionManager<T, SimplexTraitsT>::updateTableau() {
  // FIXME
  _simplexTableau.calculateRHS();
  _simplexTableau.calculateDual();
  _simplexTableau.calculateReducedCostsBasedOnDual();
  _simplexTableau.calculateSolution();
  _simplexTableau.calculateCurrentObjectiveValue();
}

template <typename T, typename SimplexTraitsT>
bool ReinversionManager<T, SimplexTraitsT>::reinverseBasisExplicit() {
  const auto basisSize = _simplexTableau._rowInfos.size();
  std::vector<int> columnIndexesMapping =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap;
  std::vector<std::vector<T>> basisColumns(basisSize);
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx)
    basisColumns[rowIdx] =
        _simplexTableau._constraintMatrixNormalForm
            ._columns[_simplexTableau.basicColumnIdx(rowIdx)];

  std::vector<std::vector<T>> newBasisMatrixInverse(basisSize);
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx) {
    newBasisMatrixInverse[rowIdx].resize(basisSize);
    newBasisMatrixInverse[rowIdx][rowIdx] = 1.0;
  }

  std::vector<bool> isUnusedColumn(basisSize, true);
  const bool isFirstVarObj = _simplexTableau._variableInfos[0]._isObjectiveVar;
  if (isFirstVarObj) {
    isUnusedColumn[0] = false;
  }
  const int firstRowIdx = ((int)isFirstVarObj);
  for (int rowIdx = firstRowIdx; rowIdx < basisSize; ++rowIdx) {
    const auto pivotColumnIdx =
        findPivotColumnMaxAbsValue(basisColumns, isUnusedColumn, rowIdx);
    if (!pivotColumnIdx.has_value()) {
      SPDLOG_WARN("BASIS MATRIX REINVERSION FAILED FOR ROW {}!", rowIdx);
      return false;
    }

    const T pivotingTermInverse{1.0 / basisColumns[*pivotColumnIdx][rowIdx]};
    const ElementaryMatrix<T> etm{basisColumns[*pivotColumnIdx],
                                  pivotingTermInverse, rowIdx};
    SimplexTraitsT::multiplyByETM(etm, newBasisMatrixInverse);

    isUnusedColumn[*pivotColumnIdx] = false;
    _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
        columnIndexesMapping[*pivotColumnIdx];

    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx])
        SimplexTraitsT::multiplyByETM(etm, basisColumns[colIdx]);
  }

  _simplexTableau._basisMatrixInverse.swap(newBasisMatrixInverse);
  return true;
}

template <typename T, typename SimplexTraitsT>
bool ReinversionManager<T, SimplexTraitsT>::reinverseBasisPFI() {
  // TODO - opt it by storing columns in vectors
  const auto basisSize = _simplexTableau._rowInfos.size();
  std::vector<int> columnIndexesMapping =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap;
  std::vector<std::vector<T>> basisColumns(basisSize);
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx)
    basisColumns[rowIdx] =
        _simplexTableau._constraintMatrixNormalForm
            ._columns[_simplexTableau.basicColumnIdx(rowIdx)];

  std::vector<bool> isUnusedColumn(basisSize, true);
  const bool isFirstVarObj = _simplexTableau._variableInfos[0]._isObjectiveVar;
  if (isFirstVarObj) {
    isUnusedColumn[0] = false;
  }
  std::vector<ElementaryMatrix<T>> newPfiEtms;
  newPfiEtms.reserve(_simplexTableau._rowInfos.size());

  const int firstRowIdx = ((int)isFirstVarObj);
  for (int rowIdx = firstRowIdx; rowIdx < basisSize; ++rowIdx) {
    const auto pivotColumnIdx =
        findPivotColumnMaxAbsValue(basisColumns, isUnusedColumn, rowIdx);
    if (!pivotColumnIdx.has_value()) {
      SPDLOG_WARN("BASIS MATRIX REINVERSION FAILED FOR ROW {}!", rowIdx);
      return false;
    }

    const T pivotingTermInverse{1.0 / basisColumns[*pivotColumnIdx][rowIdx]};
    const ElementaryMatrix<T> etm{basisColumns[*pivotColumnIdx],
                                  pivotingTermInverse, rowIdx};

    isUnusedColumn[*pivotColumnIdx] = false;
    _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
        columnIndexesMapping[*pivotColumnIdx];

    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx])
        SimplexTraitsT::multiplyByETM(etm, basisColumns[colIdx]);

    newPfiEtms.push_back(ElementaryMatrix<T>{
        std::move(basisColumns[*pivotColumnIdx]), pivotingTermInverse, rowIdx});
  }

  _simplexTableau._pfiEtms.swap(newPfiEtms);
  return true;
}

template <typename T, typename SimplexTraitsT>
bool ReinversionManager<T, SimplexTraitsT>::reinverseBasisPFISparse() {
  const auto basisSize = _simplexTableau._rowInfos.size();
  std::vector<int> columnIndexesMapping =
      _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap;
  std::vector<SparseVector<T>> basisColumns(basisSize);
  for (int rowIdx = 0; rowIdx < basisSize; ++rowIdx)
    basisColumns[rowIdx] =
        _simplexTableau._constraintMatrixSparseForm
            ._columns[_simplexTableau.basicColumnIdx(rowIdx)];

  std::vector<bool> isUnusedColumn(basisSize, true);
  const bool isFirstVarObj = _simplexTableau._variableInfos[0]._isObjectiveVar;
  if (isFirstVarObj) {
    isUnusedColumn[0] = false;
  }
  std::vector<SparseElementaryMatrix<T>> newSparsePfiEtms;
  newSparsePfiEtms.reserve(_simplexTableau._rowInfos.size());

  const int firstRowIdx = ((int)isFirstVarObj);
  for (int rowIdx = firstRowIdx; rowIdx < basisSize; ++rowIdx) {
    const auto pivotColumnIdx =
        findPivotColumnMaxAbsValue(basisColumns, isUnusedColumn, rowIdx);
    if (!pivotColumnIdx.has_value()) {
      SPDLOG_WARN("BASIS MATRIX REINVERSION FAILED FOR ROW {}!", rowIdx);
      return false;
    }

    const T pivotingTermInverse{
        1.0 / basisColumns[*pivotColumnIdx]._normalVec[rowIdx]};
    newSparsePfiEtms.push_back(SparseElementaryMatrix<T>{
        std::move(basisColumns[*pivotColumnIdx]), pivotingTermInverse, rowIdx});

    isUnusedColumn[*pivotColumnIdx] = false;
    _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx] =
        columnIndexesMapping[*pivotColumnIdx];

    for (int colIdx = 0; colIdx < basisSize; ++colIdx)
      if (isUnusedColumn[colIdx])
        SimplexTraitsT::template multiplyByETM<
            typename NumericalTraitsT::SafeNumericalAddOp>(
            newSparsePfiEtms.back(), basisColumns[colIdx]);
  }

  _simplexTableau._sparsePfiEtms.swap(newSparsePfiEtms);
  return true;
}

template <typename T, typename SimplexTraitsT>
template <typename VectorType>
std::optional<int>
ReinversionManager<T, SimplexTraitsT>::findPivotColumnFirstEligible(
    const std::vector<VectorType> &basisColumns,
    const std::vector<bool> &isUnusedColumn, const int rowIdx) const {
  const auto basisSize = _simplexTableau._rowInfos.size();
  for (int colIdx = 0; colIdx < basisSize; ++colIdx)
    if (isUnusedColumn[colIdx] &&
        NumericalTraitsT::isEligibleForPivot(basisColumns[colIdx][rowIdx]))
      return colIdx;

  return std::nullopt;
}

template <typename T, typename SimplexTraitsT>
template <typename VectorType>
std::optional<int>
ReinversionManager<T, SimplexTraitsT>::findPivotColumnMaxAbsValue(
    const std::vector<VectorType> &basisColumns,
    const std::vector<bool> &isUnusedColumn, const int rowIdx) const {
  const auto basisSize = _simplexTableau._rowInfos.size();
  std::optional<T> maxAbsPivotValue;
  std::optional<int> maxAbsPivotColIdx;
  for (int colIdx = 0; colIdx < basisSize; ++colIdx)
    if (isUnusedColumn[colIdx] &&
        NumericalTraitsT::isEligibleForPivot(basisColumns[colIdx][rowIdx])) {
      const auto currentAbsPivotValue = std::fabs(basisColumns[colIdx][rowIdx]);
      if (!maxAbsPivotColIdx.has_value() ||
          (*maxAbsPivotColIdx < currentAbsPivotValue)) {
        maxAbsPivotValue = currentAbsPivotValue;
        maxAbsPivotColIdx = colIdx;
      }
    }

  return maxAbsPivotColIdx;
}

template class ReinversionManager<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class ReinversionManager<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class ReinversionManager<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class ReinversionManager<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;