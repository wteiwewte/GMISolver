#ifndef GMISOLVER_SIMPLEXVALIDATOR_H
#define GMISOLVER_SIMPLEXVALIDATOR_H

#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/LPOptStatistics.h"

#include <unordered_set>

#include <tl/expected.hpp>

using ExpectedT = tl::expected<void, std::string>;

template <typename T, typename SimplexTraitsT> class SimplexValidator {
public:
  SimplexValidator(const SimplexTableau<T, SimplexTraitsT> &_simplexTableau,
                   const LPOptStatistics<T> &simplexLpOptStats)
      : _simplexTableau(_simplexTableau),
        _simplexLpOptStats(simplexLpOptStats) {}

  ExpectedT validatePrimalIteration(const bool isPhaseOne) const {
    return validateTableau()
        .and_then([&] { return validateBasis(); })
        .and_then([&] { return validatePrimalFeasibility(isPhaseOne); })
        .and_then(
            [&] { return validateObjectiveMonotonicity(SimplexType::PRIMAL); });
  }

  ExpectedT validateDualIteration() const {
    return validateTableau()
        .and_then([&] { return validateBasis(); })
        .and_then([&] { return validateDualFeasibility(); })
        .and_then(
            [&] { return validateObjectiveMonotonicity(SimplexType::DUAL); });
  }

  ExpectedT validateOptimality(const SimplexType simplexType,
                               const bool isPhaseOne = false) const {
    return validateTableau()
        .and_then([&] { return validateBasis(); })
        .and_then([&] { return validatePrimalFeasibility(isPhaseOne); })
        .and_then([&] { return validateDualFeasibility(); })
        .and_then([&] { return validateObjectiveMonotonicity(simplexType); });
  }

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  ExpectedT validatePrimalFeasibility(const bool isPhaseOne) const {
    return validateConstraints(isPhaseOne).and_then([&] {
      return validateVariableBounds();
    });
  }

  ExpectedT validateConstraints(const bool isPhaseOne) const {
    return validateConstraints(
               "initial_simplex", _simplexTableau._constraintMatrix,
               _simplexTableau._initialRightHandSides, _simplexTableau._x)
        .and_then([&]() -> ExpectedT {
          if (isPhaseOne)
            return {};

          return validateConstraints(
              "initial_lp", _simplexTableau._initialProgram._constraintMatrix,
              _simplexTableau._initialProgram._rightHandSides,
              _simplexTableau._x);
        });
  }

  static ExpectedT validateConstraints(const std::string &constraintsType,
                                       const Matrix<T> &constraintMatrix,
                                       const std::vector<T> &rightHandSides,
                                       const std::vector<T> &x) {
    if (constraintMatrix.empty()) {
      return {};
    }

    for (int rowIdx = 0; rowIdx < constraintMatrix.size(); ++rowIdx) {
      const int constraintMatrixWidth = constraintMatrix[0].size();
      typename SimplexTraitsT::CurrentAdder adder;
      for (int varIdx = 0; varIdx < constraintMatrixWidth; ++varIdx) {
        adder.addValue(x[varIdx] * constraintMatrix[rowIdx][varIdx]);
      }
      const auto lhs = adder.currentSum();
      const auto rhs = rightHandSides[rowIdx];
      if (!NumericalTraitsT::equal(
              lhs, rhs, NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
        return tl::unexpected{
            fmt::format("{}-th {} constraint not satisfied, diff {}", rowIdx,
                        constraintsType, std::abs(lhs - rhs))};
      }
    }
    return {};
  }

  ExpectedT validateVariableBounds() const {
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
          (_simplexTableau._x[varIdx] <
           *_simplexTableau._variableLowerBounds[varIdx] -
               NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
        return tl::unexpected{fmt::format(
            "Lower bound for var {} (is at LB {}, is at UB {}, is in basis {}) "
            "not satisfied, diff {}",
            varIdx,
            _simplexTableau._simplexBasisData
                ._isColumnAtLowerBoundBitset[varIdx],
            _simplexTableau._simplexBasisData
                ._isColumnAtUpperBoundBitset[varIdx],
            _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[varIdx],
            std::abs(_simplexTableau._x[varIdx] -
                     *_simplexTableau._variableLowerBounds[varIdx]))};
      }
      if (_simplexTableau._variableUpperBounds[varIdx].has_value() &&
          (_simplexTableau._x[varIdx] >
           *_simplexTableau._variableUpperBounds[varIdx] +
               NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
        return tl::unexpected{fmt::format(
            "Upper bound for var {} (is at LB {}, is at UB {}, is in basis {}) "
            "not satisfied, diff {}",
            varIdx,
            _simplexTableau._simplexBasisData
                ._isColumnAtLowerBoundBitset[varIdx],
            _simplexTableau._simplexBasisData
                ._isColumnAtUpperBoundBitset[varIdx],
            _simplexTableau._simplexBasisData._isBasicColumnIndexBitset[varIdx],
            std::abs(_simplexTableau._x[varIdx] -
                     *_simplexTableau._variableUpperBounds[varIdx]))};
      }
    }
    return {};
  }

  ExpectedT validateDualFeasibility() const {
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (_simplexTableau._variableInfos[varIdx]._isFixed)
        continue;

      if (_simplexTableau._variableInfos[varIdx]._isArtificial)
        continue;

      if (_simplexTableau._simplexBasisData
              ._isColumnAtLowerBoundBitset[varIdx] &&
          (_simplexTableau._reducedCosts[varIdx] <
           -NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE)) {
        return tl::unexpected{
            fmt::format("Dual feasibility not satisfied for var {} at lower "
                        "bound, reduced cost {}",
                        varIdx, _simplexTableau._reducedCosts[varIdx])};
      }

      if (_simplexTableau._simplexBasisData
              ._isColumnAtUpperBoundBitset[varIdx] &&
          (_simplexTableau._reducedCosts[varIdx] >
           NumericalTraitsT::DUAL_FEASIBILITY_TOLERANCE)) {
        return tl::unexpected{
            fmt::format("Dual feasibility not satisfied for var {} at upper "
                        "bound, reduced cost {}",
                        varIdx, _simplexTableau._reducedCosts[varIdx])};
      }
    }

    return {};
  }

  ExpectedT validateObjectiveMonotonicity(const SimplexType simplexType) const {
    const auto &consecutiveObjValues =
        _simplexLpOptStats._consecutiveObjectiveValues;
    if (consecutiveObjValues.size() <= 1)
      return {};

    const auto lastObjectiveValue =
        consecutiveObjValues[consecutiveObjValues.size() - 1];
    const auto beforeLastObjectiveValue =
        consecutiveObjValues[consecutiveObjValues.size() - 2];
    if (simplexType == SimplexType::PRIMAL) {
      if (lastObjectiveValue <
          beforeLastObjectiveValue +
              NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE)
        return {};
      else
        return tl::unexpected{fmt::format(
            "Primal obj monotonicity not satisfied - last obj value {} "
            "before last obj value {}, diff {}",
            lastObjectiveValue, beforeLastObjectiveValue,
            std::abs(lastObjectiveValue - beforeLastObjectiveValue))};
    } else {
      if (lastObjectiveValue >
          beforeLastObjectiveValue -
              NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE)
        return {};
      else
        return tl::unexpected{fmt::format(
            "Dual obj monotonicity not satisfied - last obj value {} "
            "before last obj value {}, diff {}",
            lastObjectiveValue, beforeLastObjectiveValue,
            std::abs(lastObjectiveValue - beforeLastObjectiveValue))};
    }
  }

  ExpectedT validateBasis() const {
    return validateCorrectBasisDistribution()
        .and_then([&] { return validateBoundExistence(); })
        .and_then([&] { return validateRowToColumnMapping(); });
  }

  ExpectedT validateCorrectBasisDistribution() const {
    using BitsetVecT = std::vector<boost::dynamic_bitset<>>;

    const auto &simplexBasisData = _simplexTableau._simplexBasisData;
    const auto checkIfAllDistinctBitsetPairsAreDisjoint =
        [](const BitsetVecT &bitsetVec) -> ExpectedT {
      for (int bitsetIdx1 = 0; bitsetIdx1 < bitsetVec.size(); ++bitsetIdx1) {
        for (int bitsetIdx2 = bitsetIdx1 + 1; bitsetIdx2 < bitsetVec.size();
             ++bitsetIdx2) {
          if ((bitsetVec[bitsetIdx1] & bitsetVec[bitsetIdx2]).any()) {
            return tl::unexpected{fmt::format(
                "Some variables are assigned to multiple states in basis")};
          }
        }
      }
      return {};
    };
    return checkIfAllDistinctBitsetPairsAreDisjoint(
               {simplexBasisData._isBasicColumnIndexBitset,
                simplexBasisData._isColumnAtLowerBoundBitset,
                simplexBasisData._isColumnAtUpperBoundBitset})
        .and_then([&] {
          return checkIfAllDistinctBitsetPairsAreDisjoint(
              {simplexBasisData._isColumnAtLowerBoundBitset,
               simplexBasisData._isColumnAtUpperBoundBitset,
               _simplexTableau._isVariableFreeBitset});
        })
        .and_then([&]() -> ExpectedT {
          if (!(simplexBasisData._isBasicColumnIndexBitset |
                simplexBasisData._isColumnAtLowerBoundBitset |
                simplexBasisData._isColumnAtUpperBoundBitset |
                _simplexTableau._isVariableFreeBitset)
                   .all()) {
            return tl::unexpected{fmt::format(
                "Some variables aren't assigned to any state in basis")};
          }
          return {};
        });
  }

  ExpectedT validateBoundExistence() const {
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (simplexBasisData._isColumnAtLowerBoundBitset[varIdx] &&
          !_simplexTableau._variableLowerBounds[varIdx].has_value())
        return tl::unexpected{fmt::format(
            "Var {} is assigned lower bound state but has no lower bound",
            varIdx)};

      if (simplexBasisData._isColumnAtUpperBoundBitset[varIdx] &&
          !_simplexTableau._variableUpperBounds[varIdx].has_value())
        return tl::unexpected{fmt::format(
            "Var {} is assigned upper bound state but has no upper bound",
            varIdx)};
    }
    return {};
  }

  ExpectedT validateRowToColumnMapping() const {
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;

    if (_simplexTableau._variableInfos[0]._isObjectiveVar &&
        (simplexBasisData._rowToBasisColumnIdxMap[0] != 0 ||
         !simplexBasisData._isBasicColumnIndexBitset[0]))
      return tl::unexpected{
          fmt::format("Row 0 should be mapped always to column 0 - both are "
                      "objective related")};

    std::unordered_set<int> uniqueBasicColumnIndices;
    for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
      const auto basicColumnIdx =
          simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
      if (basicColumnIdx < 0 ||
          basicColumnIdx >= _simplexTableau._variableInfos.size())
        return tl::unexpected{
            fmt::format("Basic column idx {} for row idx {} is not valid",
                        basicColumnIdx, rowIdx)};

      if (!simplexBasisData._isBasicColumnIndexBitset[basicColumnIdx])
        return tl::unexpected{
            fmt::format("Basic column idx {} has incorrect state in basis",
                        basicColumnIdx)};

      uniqueBasicColumnIndices.insert(
          simplexBasisData._rowToBasisColumnIdxMap[rowIdx]);
    }

    if (uniqueBasicColumnIndices.size() != _simplexTableau._rowInfos.size())
      return tl::unexpected{
          fmt::format("Basis row->column mapping is not unique")};

    return {};
  }

  ExpectedT validateTableau() const {
    return validateSingleVarBounds()
        .and_then([&] { return validateMatrixRepresentations(); })
        .and_then([&] { return validateRHS(); })
        .and_then([&] { return validateObjectiveRelatedThings(); })
        .and_then([&] { return validateSolutionSize(); })
        .and_then([&] { return validateBasisIntegrity(); });
  }

  ExpectedT validateSingleVarBounds() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();

    if (_simplexTableau._variableLowerBounds.size() != expectedNumberOfVars ||
        _simplexTableau._variableUpperBounds.size() != expectedNumberOfVars)
      return tl::unexpected{fmt::format(
          "Number of lower/upper bounds doesn't match number of variables")};

    if (_simplexTableau._isVariableFreeBitset.size() != expectedNumberOfVars)
      return tl::unexpected{fmt::format(
          "Size of isFreeVar bitset doesn't match number of variables")};

    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (!_simplexTableau._variableInfos[varIdx]._isFree &&
          !_simplexTableau._variableLowerBounds[varIdx].has_value()) {
        return tl::unexpected{fmt::format(
            "Non-free var {} doesn't have lower bound specified", varIdx)};
      }

      if (_simplexTableau._variableInfos[varIdx]._isFree) {
        if (_simplexTableau._variableLowerBounds[varIdx].has_value() ||
            _simplexTableau._variableUpperBounds[varIdx].has_value()) {
          return tl::unexpected{fmt::format(
              "Free var {} shouldn't have any bound specified", varIdx)};
        }
      }
    }

    return {};
  }

  ExpectedT validateMatrixRepresentations() const {
    if (_simplexTableau._simplexTableauType == SimplexTableauType::FULL) {
      return validateFullTableau();
    } else {
      const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
      const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
      if (_simplexTableau._constraintMatrix.size() != expectedNumberOfRows)
        return tl::unexpected{
            fmt::format("Number of constraints doesn't match number of rows")};

      for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
        if (_simplexTableau._constraintMatrix[rowIdx].size() !=
            expectedNumberOfVars)
          return tl::unexpected{
              fmt::format("Number of constraint row {} coefficients doesn't "
                          "match number of variables",
                          rowIdx)};
      }

      return validateMatrixReprNormal().and_then(
          [&] { return validateMatrixReprSparse(); });
    }
  }

  ExpectedT validateFullTableau() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    if (_simplexTableau._fullTableau.size() != expectedNumberOfRows)
      return tl::unexpected{fmt::format(
          "Number of full tableau constraints doesn't match number of rows")};

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._fullTableau[rowIdx].size() != expectedNumberOfVars)
        return tl::unexpected{fmt::format(
            "Number of full tableau constraint row {} coefficients doesn't "
            "match number of variables",
            rowIdx)};
    }
    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      const int mappedBasicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
      const auto &tableauRow = _simplexTableau._fullTableau[rowIdx];
      for (int varIdx = 0; varIdx < tableauRow.size(); ++varIdx) {
        if (varIdx == mappedBasicVarIdx) {
          if (!NumericalTraitsT::equal(
                  tableauRow[varIdx], 1.0,
                  NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
            return tl::unexpected{
                fmt::format("Mapped basic variable {} has value {} != 1.0 in "
                            "tableau row {}",
                            varIdx, tableauRow[varIdx], rowIdx)};
          }
        } else {
          if (_simplexTableau._simplexBasisData
                  ._isBasicColumnIndexBitset[varIdx] &&
              !NumericalTraitsT::equal(
                  tableauRow[varIdx], 0.0,
                  NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
            return tl::unexpected{fmt::format(
                "Basic variable {} has value {} != 0.0 in tableau row {} "
                "mapped to variable {}",
                varIdx, tableauRow[varIdx], rowIdx, mappedBasicVarIdx)};
          }
        }
      }
    }

    return {};
  }

  ExpectedT validateMatrixReprNormal() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();

    if (_simplexTableau._constraintMatrixNormalForm._rows.size() !=
        expectedNumberOfRows)
      return tl::unexpected{fmt::format(
          "[Normal form] Number of constraints doesn't match number of rows")};

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._constraintMatrixNormalForm._rows[rowIdx].size() !=
          expectedNumberOfVars)
        return tl::unexpected{
            fmt::format("[Normal form] Number of constraint row {} "
                        "coefficients doesn't match number of variables",
                        rowIdx)};
    }

    if (_simplexTableau._constraintMatrixNormalForm._columns.size() !=
        expectedNumberOfVars)
      return tl::unexpected{
          fmt::format("[Normal form] Number of constraints columns doesn't "
                      "match number of variables")};

    for (int colIdx = 0; colIdx < expectedNumberOfVars; ++colIdx) {
      if (_simplexTableau._constraintMatrixNormalForm._columns[colIdx].size() !=
          expectedNumberOfRows)
        return tl::unexpected{
            fmt::format("[Normal form] Number of constraint column {} "
                        "coefficients doesn't match number of variables",
                        colIdx)};
    }

    return {};
  }

  ExpectedT validateMatrixReprSparse() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();

    if (_simplexTableau._constraintMatrixSparseForm._rows.size() !=
        expectedNumberOfRows)
      return tl::unexpected{fmt::format(
          "[Sparse form] Number of constraints doesn't match number of rows")};

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._constraintMatrixSparseForm._rows[rowIdx]
              ._normalVec.size() != expectedNumberOfVars)
        return tl::unexpected{
            fmt::format("[Sparse form] Number of constraint row {} "
                        "coefficients doesn't match number of variables",
                        rowIdx)};
    }

    if (_simplexTableau._constraintMatrixSparseForm._columns.size() !=
        expectedNumberOfVars)
      return tl::unexpected{
          fmt::format("[Sparse form] Number of constraints columns doesn't "
                      "match number of variables")};

    for (int colIdx = 0; colIdx < expectedNumberOfVars; ++colIdx) {
      if (_simplexTableau._constraintMatrixSparseForm._columns[colIdx]
              ._normalVec.size() != expectedNumberOfRows)
        return tl::unexpected{
            fmt::format("[Sparse form] Number of constraint column {} "
                        "coefficients doesn't match number of variables",
                        colIdx)};
    }

    return {};
  }

  ExpectedT validateRHS() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    if (_simplexTableau._rightHandSides.size() != expectedNumberOfRows)
      return tl::unexpected{
          fmt::format("RHS size doesn't match number of rows")};

    if (_simplexTableau._simplexTableauType != SimplexTableauType::FULL) {
      if (_simplexTableau._initialRightHandSides.size() != expectedNumberOfRows)
        return tl::unexpected{
            fmt::format("Initial RHS size doesn't match number of rows")};
    }

    return {};
  }

  ExpectedT validateObjectiveRelatedThings() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    if (_simplexTableau._objectiveRow.size() == expectedNumberOfVars &&
        _simplexTableau._reducedCosts.size() == expectedNumberOfVars)
      return {};

    return tl::unexpected{fmt::format("Objective row or reduced costs sizes "
                                      "don't match number of variables")};
  }

  ExpectedT validateSolutionSize() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    if (_simplexTableau._x.size() == expectedNumberOfVars)
      return {};

    return tl::unexpected{fmt::format(
        "Solution vector size {} doesn't match number of variables {}",
        _simplexTableau._x.size(), expectedNumberOfVars)};
  }

  ExpectedT validateBasisIntegrity() const {
    return validateBasisDataSizes().and_then(
        [&] { return validateBasisReprSizes(); });
  }

  ExpectedT validateBasisDataSizes() const {
    const auto expectedNumberOfVars = _simplexTableau._variableInfos.size();
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;

    if (simplexBasisData._rowToBasisColumnIdxMap.size() ==
            expectedNumberOfRows &&
        simplexBasisData._isBasicColumnIndexBitset.size() ==
            expectedNumberOfVars &&
        simplexBasisData._isColumnAtLowerBoundBitset.size() ==
            expectedNumberOfVars &&
        simplexBasisData._isColumnAtUpperBoundBitset.size() ==
            expectedNumberOfVars)
      return {};

    return tl::unexpected{fmt::format(
        "Basis state vector sizes don't match number of rows/variables")};
  }

  ExpectedT validateBasisReprSizes() const {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue) {
      if (_simplexTableau._simplexTableauType !=
          SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE)
        return tl::unexpected{
            fmt::format("Sparse representation is valid only for PFI mode")};

      return validatePFISparse();
    }

    switch (_simplexTableau._simplexTableauType) {
    case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
      return validatePFINormal();
    }
    case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
      return validateBasisMatrixInverse();
    default:
      return {};
    }
  }

  ExpectedT validatePFINormal() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    if (std::all_of(_simplexTableau._pfiEtms.begin(),
                    _simplexTableau._pfiEtms.end(),
                    [&](const ElementaryMatrix<T> &etm) {
                      return etm._vec.size() == expectedNumberOfRows;
                    }))
      return {};

    return tl::unexpected{fmt::format(
        "[Normal form] Some of PFI etms sizes don't match number of rows")};
  }

  ExpectedT validatePFISparse() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    if (std::all_of(_simplexTableau._sparsePfiEtms.begin(),
                    _simplexTableau._sparsePfiEtms.end(),
                    [&](const SparseElementaryMatrix<T> &etm) {
                      return etm._sparseVec._normalVec.size() ==
                             expectedNumberOfRows;
                    }))
      return {};

    return tl::unexpected{fmt::format(
        "[Sparse form] Some of PFI etms sizes don't match number of rows")};
  }

  ExpectedT validateBasisMatrixInverse() const {
    const auto expectedNumberOfRows = _simplexTableau._rowInfos.size();
    if (_simplexTableau._basisMatrixInverse.size() != expectedNumberOfRows)
      return tl::unexpected{fmt::format(
          "Basis matrix inverse size doesn't match number of rows")};

    for (int rowIdx = 0; rowIdx < expectedNumberOfRows; ++rowIdx) {
      if (_simplexTableau._basisMatrixInverse[rowIdx].size() !=
          expectedNumberOfRows)
        return tl::unexpected{fmt::format(
            "Basis matrix inverse row {} size doesn't match number of rows",
            rowIdx)};
    }

    return {};
  }

  const SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const LPOptStatistics<T> &_simplexLpOptStats;
};

#endif // GMISOLVER_SIMPLEXVALIDATOR_H
