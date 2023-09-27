#ifndef GMISOLVER_SIMPLEXVALIDATOR_H
#define GMISOLVER_SIMPLEXVALIDATOR_H

#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/LPOptStatistics.h"

#include <unordered_set>

template <typename T, typename SimplexTraitsT> class SimplexValidator {
public:
  SimplexValidator(const SimplexTableau<T, SimplexTraitsT> &_simplexTableau,
                   const LPOptStatistics<T> &simplexLpOptStats)
      : _simplexTableau(_simplexTableau),
        _simplexLpOptStats(simplexLpOptStats) {}

  bool validatePrimalIteration() const {
    return validateTableau() && validateBasis() &&
           validatePrimalFeasibility() &&
           validateObjectiveMonotonicity(SimplexType::PRIMAL);
  }

  bool validateDualIteration() const {
    return validateTableau() && validateBasis() && validateDualFeasibility() &&
           validateObjectiveMonotonicity(SimplexType::DUAL);
  }

  bool validateOptimality(const SimplexType simplexType) const {
    return validateTableau() && validateBasis() &&
           validatePrimalFeasibility() && validateDualFeasibility() &&
           validateObjectiveMonotonicity(simplexType);
  }

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  bool validatePrimalFeasibility() const {
    return validateConstraints() && validateVariableBounds();
  }

  bool validateConstraints() const {
    for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
      typename SimplexTraitsT::CurrentAdder adder;
      for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
           ++varIdx) {
        adder.addValue(_simplexTableau._x[varIdx] *
                       _simplexTableau._constraintMatrix[rowIdx][varIdx]);
      }
      const auto lhs = adder.currentSum();
      const auto rhs = _simplexTableau._initialRightHandSides[rowIdx];
      if (!NumericalTraitsT::equal(
              lhs, rhs, NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
        return false;
      }
    }
    return true;
  }

  bool validateVariableBounds() const {
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
          (_simplexTableau._x[varIdx] <
           *_simplexTableau._variableLowerBounds[varIdx] -
               NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
        return false;
      }
      if (_simplexTableau._variableUpperBounds[varIdx].has_value() &&
          (_simplexTableau._x[varIdx] >
           *_simplexTableau._variableUpperBounds[varIdx] +
               NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
        return false;
      }
    }
    return true;
  }

  bool validateDualFeasibility() const {
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (_simplexTableau._variableInfos[varIdx]._isFixed)
        continue;

      if (_simplexTableau._simplexBasisData
              ._isColumnAtLowerBoundBitset[varIdx] &&
          (_simplexTableau._reducedCosts[varIdx] <
           -NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE)) {
        return false;
      }

      if (_simplexTableau._simplexBasisData
              ._isColumnAtUpperBoundBitset[varIdx] &&
          (_simplexTableau._reducedCosts[varIdx] >
           NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE)) {
        return false;
      }
    }

    return true;
  }

  bool validateObjectiveMonotonicity(const SimplexType simplexType) const {
    const auto &consecutiveObjValues =
        _simplexLpOptStats._consecutiveObjectiveValues;
    if (consecutiveObjValues.size() <= 1)
      return true;

    const auto lastObjectiveValue =
        consecutiveObjValues[consecutiveObjValues.size() - 1];
    const auto beforeLastObjectiveValue =
        consecutiveObjValues[consecutiveObjValues.size() - 2];
    if (simplexType == SimplexType::PRIMAL)
      return lastObjectiveValue <
             beforeLastObjectiveValue +
                 NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE;
    else
      return lastObjectiveValue >
             beforeLastObjectiveValue -
                 NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE;
  }

  bool validateBasis() const {
    return validateCorrectBasisDistribution() && validateBoundExistence() &&
           validateRowToColumnMapping();
  }

  bool validateCorrectBasisDistribution() const {
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;
    if ((simplexBasisData._isBasicColumnIndexBitset &
         simplexBasisData._isColumnAtLowerBoundBitset &
         simplexBasisData._isColumnAtUpperBoundBitset)
            .any()) {
      return false;
    }
    if (!(simplexBasisData._isBasicColumnIndexBitset |
          simplexBasisData._isColumnAtLowerBoundBitset |
          simplexBasisData._isColumnAtUpperBoundBitset)
             .all()) {
      return false;
    }
    return true;
  }

  bool validateBoundExistence() const {
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (simplexBasisData._isColumnAtLowerBoundBitset[varIdx] &&
          !_simplexTableau._variableLowerBounds[varIdx].has_value())
        return false;

      if (simplexBasisData._isColumnAtUpperBoundBitset[varIdx] &&
          !_simplexTableau._variableUpperBounds[varIdx].has_value())
        return false;
    }
    return true;
  }

  bool validateRowToColumnMapping() const {
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;
    std::unordered_set<int> uniqueBasicColumnIndices;
    for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
      const auto basicColumnIdx =
          simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
      if (basicColumnIdx < 0 ||
          basicColumnIdx >= _simplexTableau._variableInfos.size())
        return false;

      uniqueBasicColumnIndices.insert(
          simplexBasisData._rowToBasisColumnIdxMap[rowIdx]);
    }

    if (uniqueBasicColumnIndices.size() != _simplexTableau._rowInfos.size())
      return false;

    return true;
  }

  bool validateTableau() const {
    return validateSingleVarBounds() && validateMatrixRepresentations() &&
           validateRHS() && validateObjectiveRelatedThings() &&
           validateSolutionSize() && validateBasisIntegrity();
  }

  bool validateSingleVarBounds() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    if (_simplexTableau._variableLowerBounds.size() != expectedNumberOfVars ||
        _simplexTableau._variableUpperBounds.size() != expectedNumberOfVars)
      return false;

    return std::all_of(
        _simplexTableau._variableLowerBounds.begin(),
        _simplexTableau._variableLowerBounds.end(),
        [](const std::optional<T> &lb) { return lb.has_value(); });
  }

  bool validateMatrixRepresentations() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();

    if (_simplexTableau._constraintMatrix.size() != expectedNumberOfRows)
      return false;

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._constraintMatrix[rowIdx].size() !=
          expectedNumberOfVars)
        return false;
    }

    return validateMatrixReprNormal() && validateMatrixReprSparse();
  }

  bool validateMatrixReprNormal() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();

    if (_simplexTableau._constraintMatrixNormalForm._rows.size() !=
        expectedNumberOfRows)
      return false;

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._constraintMatrixNormalForm._rows[rowIdx].size() !=
          expectedNumberOfVars)
        return false;
    }

    if (_simplexTableau._constraintMatrixNormalForm._columns.size() !=
        expectedNumberOfVars)
      return false;

    for (int colIdx = 0; colIdx < expectedNumberOfVars; ++colIdx) {
      if (_simplexTableau._constraintMatrixNormalForm._columns[colIdx].size() !=
          expectedNumberOfRows)
        return false;
    }

    return true;
  }

  bool validateMatrixReprSparse() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();

    if (_simplexTableau._constraintMatrixSparseForm._rows.size() !=
        expectedNumberOfRows)
      return false;

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._constraintMatrixSparseForm._rows[rowIdx]
              ._normalVec.size() != expectedNumberOfVars)
        return false;
    }

    if (_simplexTableau._constraintMatrixSparseForm._columns.size() !=
        expectedNumberOfVars)
      return false;

    for (int colIdx = 0; colIdx < expectedNumberOfVars; ++colIdx) {
      if (_simplexTableau._constraintMatrixSparseForm._columns[colIdx]
              ._normalVec.size() != expectedNumberOfRows)
        return false;
    }

    return true;
  }

  bool validateRHS() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    return _simplexTableau._rightHandSides.size() == expectedNumberOfRows &&
           _simplexTableau._initialRightHandSides.size() ==
               expectedNumberOfRows;
  }

  bool validateObjectiveRelatedThings() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    return _simplexTableau._objectiveRow.size() == expectedNumberOfVars &&
           _simplexTableau._reducedCosts.size() == expectedNumberOfVars;
  }

  bool validateSolutionSize() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    return _simplexTableau._x.size() == expectedNumberOfVars;
  }

  bool validateBasisIntegrity() const {
    return validateBasisDataSizes() && validateBasisReprSizes();
  }

  bool validateBasisDataSizes() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;

    return simplexBasisData._rowToBasisColumnIdxMap.size() ==
               expectedNumberOfRows &&
           simplexBasisData._isBasicColumnIndexBitset.size() ==
               expectedNumberOfVars &&
           simplexBasisData._isColumnAtLowerBoundBitset.size() ==
               expectedNumberOfVars &&
           simplexBasisData._isColumnAtUpperBoundBitset.size() ==
               expectedNumberOfVars;
  }

  bool validateBasisReprSizes() const {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue) {
      if (!_simplexTableau._useProductFormOfInverse)
        return false;

      return validatePFISparse();
    }

    return _simplexTableau._useProductFormOfInverse
               ? validatePFINormal()
               : validateBasisMatrixInverse();
  }

  bool validatePFINormal() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    return std::all_of(_simplexTableau._pfiEtms.begin(),
                       _simplexTableau._pfiEtms.end(),
                       [&](const ElementaryMatrix<T> &etm) {
                         return etm._vec.size() == expectedNumberOfRows;
                       });
  }

  bool validatePFISparse() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    return std::all_of(_simplexTableau._sparsePfiEtms.begin(),
                       _simplexTableau._sparsePfiEtms.end(),
                       [&](const SparseElementaryMatrix<T> &etm) {
                         return etm._sparseVec._normalVec.size() ==
                                expectedNumberOfRows;
                       });
  }

  bool validateBasisMatrixInverse() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    if (_simplexTableau._basisMatrixInverse.size() != expectedNumberOfRows)
      return false;

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._basisMatrixInverse[rowIdx].size() !=
          expectedNumberOfRows)
        return false;
    }

    return true;
  }

  const SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const LPOptStatistics<T> &_simplexLpOptStats;
};

#endif // GMISOLVER_SIMPLEXVALIDATOR_H
