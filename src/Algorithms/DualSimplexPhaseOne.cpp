#include "src/Algorithms/DualSimplexPhaseOne.h"

#include "src/Algorithms/DualSimplex.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexTableauResizer.h"
#include "src/DataModel/LinearProgram.h"
#include "src/Util/SpdlogHeader.h"
#include "src/Util/Time.h"

template <typename T, typename SimplexTraitsT>
DualSimplexPhaseOne<T, SimplexTraitsT>::DualSimplexPhaseOne(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const DualSimplexRowPivotRule dualSimplexRowPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager),
      _dualSimplexRowPivotRule(dualSimplexRowPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption) {

  _simplexTableau.addArtificialVariables(SimplexType::DUAL);
  _simplexTableau.initMatrixRepresentations();

  _simplexTableau.init(SimplexType::DUAL);
  if (_simplexTableau._simplexTableauType == SimplexTableauType::FULL) {
    _simplexTableau.initBasisMatrixInverse();
  }
  markBoxedVariablesAsNotEligible();
  changeBoundsRHSAndObj();

  _simplexTableau.setObjective(auxiliaryObjective());
  initBoundsForBoxedVars();

  _simplexTableau.calculateRHS();
  _simplexTableau.calculateSolution();
  _simplexTableau.calculateCurrentObjectiveValue();
  SPDLOG_TRACE("Simplex tableau with artificial variables");
  SPDLOG_TRACE(toStringLpSolveFormat());
}

template <typename T, typename SimplexTraitsT>
std::string DualSimplexPhaseOne<T, SimplexTraitsT>::type() const {
  return fmt::format(
      "DUAL SIMPLEX PHASE ONE ({}, {})",
      simplexTableauTypeToStr(_simplexTableau._simplexTableauType),
      matrixRepresentationTypeToStr(SimplexTraitsT::matrixRepresentationType));
}

template <typename T, typename SimplexTraitsT>
void DualSimplexPhaseOne<T, SimplexTraitsT>::initBoundsForBoxedVars() {
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (!_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[varIdx] &&
        _simplexTableau._variableLowerBounds[varIdx].has_value() &&
        _simplexTableau._variableUpperBounds[varIdx].has_value()) {
      if (_simplexTableau._reducedCosts[varIdx] < 0.0) {
        if (!_simplexTableau._variableUpperBounds[varIdx].has_value()) {
          SPDLOG_ERROR("VAR (IDX {}, LABEL {}) HAS NEGATIVE REDUCED COST {}, "
                       "BUT DOESN'T HAVE "
                       "UPPER BOUND",
                       varIdx, _simplexTableau._variableInfos[varIdx]._label,
                       _simplexTableau._reducedCosts[varIdx]);
          continue;
        }
        _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[varIdx] =
            false;
        _simplexTableau._simplexBasisData._isColumnAtUpperBoundBitset[varIdx] =
            true;
      } else {
        if (!_simplexTableau._variableLowerBounds[varIdx].has_value()) {
          SPDLOG_ERROR("VAR (IDX {}, LABEL {}) HAS NONNEGATIVE REDUCED COST "
                       "{}, BUT DOESN'T HAVE "
                       "LOWER BOUND",
                       varIdx, _simplexTableau._variableInfos[varIdx]._label,
                       _simplexTableau._reducedCosts[varIdx]);
          continue;
        }
        _simplexTableau._simplexBasisData._isColumnAtLowerBoundBitset[varIdx] =
            true;
        _simplexTableau._simplexBasisData._isColumnAtUpperBoundBitset[varIdx] =
            false;
      }
    }
  }
}

template <typename T, typename SimplexTraitsT>
void DualSimplexPhaseOne<T, SimplexTraitsT>::markBoxedVariablesAsNotEligible() {
  auto &isColumnEligibleBitset =
      _simplexTableau._simplexBasisData._isColumnEligibleBitset;
  isColumnEligibleBitset = boost::dynamic_bitset<>{};
  isColumnEligibleBitset.value().resize(_simplexTableau._variableInfos.size(),
                                        true);
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (_simplexTableau._variableInfos[varIdx]._isArtificial)
      continue;

    if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
        _simplexTableau._variableUpperBounds[varIdx].has_value()) {
      isColumnEligibleBitset.value()[varIdx] = false;
    }
  }
}

template <typename T, typename SimplexTraitsT>
void DualSimplexPhaseOne<T,
                         SimplexTraitsT>::unmarkBoxedVariablesAsNotEligible() {
  _simplexTableau._simplexBasisData._isColumnEligibleBitset = std::nullopt;
}

template <typename T, typename SimplexTraitsT>
void DualSimplexPhaseOne<T, SimplexTraitsT>::changeBoundsRHSAndObj() {
  _originalVariableInfos = _simplexTableau._variableInfos;
  _originalVariableLowerBounds = _simplexTableau._variableLowerBounds;
  _originalVariableUpperBounds = _simplexTableau._variableUpperBounds;
  _originalInitialRightHandSides = _simplexTableau._initialRightHandSides;

  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
        _simplexTableau._variableUpperBounds[varIdx].has_value()) {
      _simplexTableau._variableLowerBounds[varIdx] = 0;
      _simplexTableau._variableUpperBounds[varIdx] = 0;
      _simplexTableau._variableInfos[varIdx]._isFixed = true;
    } else {
      if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
          !_simplexTableau._variableUpperBounds[varIdx].has_value()) {
        _simplexTableau._variableLowerBounds[varIdx] = 0;
        _simplexTableau._variableUpperBounds[varIdx] = 1;
      } else if (!_simplexTableau._variableLowerBounds[varIdx].has_value() &&
                 _simplexTableau._variableUpperBounds[varIdx].has_value()) {
        _simplexTableau._variableLowerBounds[varIdx] = -1;
        _simplexTableau._variableUpperBounds[varIdx] = 0;
      } else if (!_simplexTableau._variableLowerBounds[varIdx].has_value() &&
                 !_simplexTableau._variableUpperBounds[varIdx].has_value()) {
        _simplexTableau._variableLowerBounds[varIdx] = -1;
        _simplexTableau._variableUpperBounds[varIdx] = 1;
      }
    }
    _simplexTableau._variableInfos[varIdx]._isFree = false;
  }
  if (_simplexTableau._variableInfos[0]._isObjectiveVar) {
    _simplexTableau._variableInfos[0]._isObjectiveVar = false;
  }
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    _simplexTableau._initialRightHandSides[rowIdx] = 0;
  }
}

template <typename T, typename SimplexTraitsT>
void DualSimplexPhaseOne<T, SimplexTraitsT>::restoreBoundsRHSAndObj() {
  _simplexTableau._variableLowerBounds = _originalVariableLowerBounds;
  _simplexTableau._variableUpperBounds = _originalVariableUpperBounds;
  _simplexTableau._initialRightHandSides = _originalInitialRightHandSides;
  _simplexTableau._variableInfos = _originalVariableInfos;

  _simplexTableau.setObjective(_simplexTableau.originalObjective());
  initBoundsForBoxedVars();
  recalculateRHS();
  _simplexTableau.calculateSolution();
  _simplexTableau.calculateCurrentObjectiveValue();
  if (_simplexTableau._simplexTableauType == SimplexTableauType::FULL) {
    _simplexTableau._basisMatrixInverse.clear();
  }
}

template <typename T, typename SimplexTraitsT>
void DualSimplexPhaseOne<T, SimplexTraitsT>::recalculateRHS() {
  switch (_simplexTableau._simplexTableauType) {
  case SimplexTableauType::FULL:
    return _simplexTableau.calculateRHSBasisInverse();
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE: {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue)
      return _simplexTableau.calculateRHSPFISparse();

    return _simplexTableau.calculateRHSPFI();
  }
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return _simplexTableau.calculateRHSBasisInverse();
  }
}

template <typename T, typename SimplexTraitsT>
std::vector<T>
DualSimplexPhaseOne<T, SimplexTraitsT>::auxiliaryObjective() const {
  std::vector<T> auxObjective = _simplexTableau.originalObjective();
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
        _simplexTableau._variableUpperBounds[varIdx].has_value()) {
      auxObjective[varIdx] = 0;
    }
  }

  return auxObjective;
}

template <typename T, typename SimplexTraitsT>
DualSimplex<T, SimplexTraitsT>
DualSimplexPhaseOne<T, SimplexTraitsT>::dualSimplex() const {
  return DualSimplex<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _dualSimplexRowPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T> DualSimplexPhaseOne<T, SimplexTraitsT>::run() {
  auto phaseOneLPStats = dualSimplex().run(
      "DUAL_PHASE_ONE", PrintSimplexOptSummary::YES, DualPhase::ONE);

  if (phaseOneLPStats._optResult !=
      LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
    SPDLOG_WARN("MIN SUM OF INFEASIBILITIES PROGRAM RETURNED {} RESULT",
                lpOptimizationResultToStr(phaseOneLPStats._optResult));
    phaseOneLPStats._phaseOneSucceeded = false;
    phaseOneLPStats._optResult = LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED;
    return phaseOneLPStats;
  }

  if (_simplexTableau._objectiveValue <
      -NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE) {
    SPDLOG_WARN(
        "MIN SUM OF INFEASIBILITIES PROGRAM HAS OPTIMUM OF {} LESS THAN 0 - "
        "INITIAL PROGRAM IS DUAL INFEASIBLE",
        _simplexTableau._objectiveValue);
    phaseOneLPStats._phaseOneSucceeded = false;
    phaseOneLPStats._optResult = LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED;
    return phaseOneLPStats;
  }

  phaseOneLPStats._phaseOneSucceeded = true;
  unmarkBoxedVariablesAsNotEligible();
  restoreBoundsRHSAndObj();
  return phaseOneLPStats;
}

template class DualSimplexPhaseOne<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexPhaseOne<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class DualSimplexPhaseOne<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexPhaseOne<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
