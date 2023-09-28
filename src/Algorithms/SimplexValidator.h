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

  ExpectedT validatePrimalIteration() const {
    return validateTableau()
        .and_then([&] { return validateBasis(); })
        .and_then([&] { return validatePrimalFeasibility(); })
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

  ExpectedT validateOptimality(const SimplexType simplexType) const {
    return validateTableau()
        .and_then([&] { return validateBasis(); })
        .and_then([&] { return validatePrimalFeasibility(); })
        .and_then([&] { return validateDualFeasibility(); })
        .and_then([&] { return validateObjectiveMonotonicity(simplexType); });
  }

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  ExpectedT validatePrimalFeasibility() const {
    return validateConstraints().and_then(
        [&] { return validateVariableBounds(); });
  }

  ExpectedT validateConstraints() const {
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
        return tl::unexpected{
            fmt::format("{}-th constraint not satisfied, diff {}", rowIdx,
                        std::abs(lhs - rhs))};
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
            "Lower bound for var {} not satisfied, diff {}", varIdx,
            std::abs(_simplexTableau._x[varIdx] -
                     *_simplexTableau._variableLowerBounds[varIdx]))};
      }
      if (_simplexTableau._variableUpperBounds[varIdx].has_value() &&
          (_simplexTableau._x[varIdx] >
           *_simplexTableau._variableUpperBounds[varIdx] +
               NumericalTraitsT::PRIMAL_FEASIBILITY_TOLERANCE)) {
        return tl::unexpected{fmt::format(
            "Upper bound for var {} not satisfied, diff {}", varIdx,
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
    const auto &simplexBasisData = _simplexTableau._simplexBasisData;
    if ((simplexBasisData._isBasicColumnIndexBitset &
         simplexBasisData._isColumnAtLowerBoundBitset &
         simplexBasisData._isColumnAtUpperBoundBitset)
            .any()) {
      return tl::unexpected{fmt::format(
          "Some variables are assigned to multiple states in basis")};
    }
    if (!(simplexBasisData._isBasicColumnIndexBitset |
          simplexBasisData._isColumnAtLowerBoundBitset |
          simplexBasisData._isColumnAtUpperBoundBitset)
             .all()) {
      return tl::unexpected{
          fmt::format("Some variables aren't assigned to any state in basis")};
    }
    return {};
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

    if (!std::all_of(_simplexTableau._variableLowerBounds.begin(),
                     _simplexTableau._variableLowerBounds.end(),
                     [](const std::optional<T> &lb) { return lb.has_value(); }))
      return tl::unexpected{
          fmt::format("Some of variables don't have lower bound specified")};

    return {};
  }

  ExpectedT validateMatrixRepresentations() const {
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
    if (_simplexTableau._rightHandSides.size() == expectedNumberOfRows &&
        _simplexTableau._initialRightHandSides.size() == expectedNumberOfRows)
      return {};

    return tl::unexpected{
        fmt::format("RHS or initial RHS sizes don't match number of rows")};
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

    return tl::unexpected{
        fmt::format("Solution vector size doesn't match number of variables")};
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
      if (!_simplexTableau._useProductFormOfInverse)
        return tl::unexpected{
            fmt::format("Sparse representation is valid only for PFI mode")};

      return validatePFISparse();
    }

    return _simplexTableau._useProductFormOfInverse
               ? validatePFINormal()
               : validateBasisMatrixInverse();
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
