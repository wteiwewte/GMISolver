#include "src/Algorithms/PhaseOneUtilities.h"

#include "src/Algorithms/SimplexTableau.h"

template <typename T, typename SimplexTraitsT>
PhaseOneUtilities<T, SimplexTraitsT>::PhaseOneUtilities(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau)
    : _simplexTableau(simplexTableau) {}

template <typename T, typename SimplexTraitsT>
void PhaseOneUtilities<T, SimplexTraitsT>::addArtificialVariables(
    const SimplexType simplexType) {
  const int variableCountAtTheStart = _simplexTableau._variableInfos.size();
  const bool shouldMapZerothVarToZerothRow =
      shouldMapZerothVarToZerothRowImpl(simplexType);
  const int artificialVarCount =
      _simplexTableau._rowInfos.size() - ((int)shouldMapZerothVarToZerothRow);
  const int newVariableCount = variableCountAtTheStart + artificialVarCount;

  SPDLOG_INFO("ADDING {} ARTIFICIAL VARIABLES", artificialVarCount);

  const auto newArtificialLabel = [&](const auto varIdx) {
    return "A" + std::to_string(varIdx);
  };

  _simplexTableau._constraintMatrix[0].resize(newVariableCount);
  const int firstRowIdx = ((int)shouldMapZerothVarToZerothRow);
  for (int rowIdx = firstRowIdx; rowIdx < _simplexTableau._rowInfos.size();
       ++rowIdx) {
    _simplexTableau._constraintMatrix[rowIdx].resize(newVariableCount);

    const auto newVariableIdx = variableCountAtTheStart + rowIdx - firstRowIdx;
    _simplexTableau._constraintMatrix[rowIdx][newVariableIdx] = 1;
    const auto newArtificialLabelStr = newArtificialLabel(newVariableIdx);
    if (simplexType == SimplexType::PRIMAL) {
      _simplexTableau._variableLowerBounds.push_back(0.0);
      _simplexTableau._variableUpperBounds.push_back(std::nullopt);
      _simplexTableau._variableInfos.push_back(
          VariableInfo{._label = newArtificialLabelStr,
                       ._type = VariableType::CONTINUOUS,
                       ._isArtificial = true});
    } else {
      _simplexTableau._variableLowerBounds.push_back(0.0);
      _simplexTableau._variableUpperBounds.push_back(0.0);
      _simplexTableau._variableInfos.push_back(
          VariableInfo{._label = newArtificialLabelStr,
                       ._type = VariableType::INTEGER,
                       ._isArtificial = true,
                       ._isFixed = true});
    }
  }
}

template <typename T, typename SimplexTraitsT>
std::vector<T>
PhaseOneUtilities<T, SimplexTraitsT>::artificialObjective() const {
  std::vector<T> result(_simplexTableau._constraintMatrix.front().size());

  for (int variableIdx = 0; variableIdx < _simplexTableau._variableInfos.size();
       ++variableIdx)
    if (_simplexTableau._variableInfos[variableIdx]._isArtificial)
      result[variableIdx] = 1;

  return result;
}

template <typename T, typename SimplexTraitsT>
std::optional<SimplexBasisData>
PhaseOneUtilities<T, SimplexTraitsT>::createBasisFromArtificialVars(
    const SimplexType simplexType) const {
  std::optional<int> firstArtificialIdx;
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size(); ++varIdx)
    if (_simplexTableau._variableInfos[varIdx]._isArtificial) {
      firstArtificialIdx = varIdx;
      break;
    }

  if (!firstArtificialIdx.has_value()) {
    SPDLOG_ERROR("No artificial variable found");
    return std::nullopt;
  }

  SimplexBasisData result;
  result._isBasicColumnIndexBitset.resize(
      _simplexTableau._variableInfos.size());
  result._isColumnAtLowerBoundBitset.resize(
      _simplexTableau._variableInfos.size(), false);
  result._isColumnAtUpperBoundBitset.resize(
      _simplexTableau._variableInfos.size(), false);
  result._rowToBasisColumnIdxMap.resize(_simplexTableau._rowInfos.size());

  const bool shouldMapZerothVarToZerothRow =
      shouldMapZerothVarToZerothRowImpl(simplexType);
  if (shouldMapZerothVarToZerothRow) {
    result._rowToBasisColumnIdxMap[0] = 0;
    result._isBasicColumnIndexBitset[0] = true;
  }

  const int firstVarIdx = ((int)shouldMapZerothVarToZerothRow);
  for (int varIdx = firstVarIdx; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    result._isColumnAtLowerBoundBitset[varIdx] =
        !_simplexTableau._variableInfos[varIdx]._isFree;
  }

  const int firstRowIdx = ((int)shouldMapZerothVarToZerothRow);
  for (int rowIdx = firstRowIdx; rowIdx < _simplexTableau._rowInfos.size();
       ++rowIdx) {
    const auto basicColumnIdx = *firstArtificialIdx + rowIdx - firstRowIdx;
    result._rowToBasisColumnIdxMap[rowIdx] = basicColumnIdx;
    result._isBasicColumnIndexBitset[basicColumnIdx] = true;
    result._isColumnAtLowerBoundBitset[basicColumnIdx] = false;
  }

  return result;
}

template <typename T, typename SimplexTraitsT>
bool PhaseOneUtilities<T, SimplexTraitsT>::shouldMapZerothVarToZerothRowImpl(
    const SimplexType simplexType) const {
  return simplexType == SimplexType::PRIMAL &&
         _simplexTableau._variableInfos[0]._isObjectiveVar;
}

template class PhaseOneUtilities<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class PhaseOneUtilities<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class PhaseOneUtilities<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class PhaseOneUtilities<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
