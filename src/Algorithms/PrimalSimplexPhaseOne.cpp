#include "src/Algorithms/PrimalSimplexPhaseOne.h"

#include "src/Algorithms/PrimalSimplex.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexTableauResizer.h"
#include "src/DataModel/LinearProgram.h"
#include "src/Util/SpdlogHeader.h"
#include "src/Util/Time.h"

template <typename T, typename SimplexTraitsT>
PrimalSimplexPhaseOne<T, SimplexTraitsT>::PrimalSimplexPhaseOne(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption) {
  makeRightHandSidesNonNegative();
  _simplexTableau.addArtificialVariables(SimplexType::PRIMAL);
  _simplexTableau.initMatrixRepresentations();
  _simplexTableau.init(SimplexType::PRIMAL);
  if (_simplexTableau._simplexTableauType == SimplexTableauType::FULL) {
    _simplexTableau.initBasisMatrixInverse();
  }

  _simplexTableau.setObjective(_simplexTableau.artificialObjective());

  _simplexTableau.calculateRHS();
  _simplexTableau.calculateSolution();
  _simplexTableau.calculateCurrentObjectiveValue();
  SPDLOG_TRACE("Simplex tableau with artificial variables");
  SPDLOG_TRACE(toString());
  SPDLOG_TRACE(toStringLpSolveFormat());
}

template <typename T, typename SimplexTraitsT>
std::string PrimalSimplexPhaseOne<T, SimplexTraitsT>::type() const {
  return fmt::format(
      "PRIMAL SIMPLEX PHASE ONE ({}, {})",
      simplexTableauTypeToStr(_simplexTableau._simplexTableauType),
      matrixRepresentationTypeToStr(SimplexTraitsT::matrixRepresentationType));
}

template <typename T, typename SimplexTraitsT>
PrimalSimplex<T, SimplexTraitsT>
PrimalSimplexPhaseOne<T, SimplexTraitsT>::primalSimplex() const {
  return PrimalSimplex<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _primalSimplexColumnPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T> PrimalSimplexPhaseOne<T, SimplexTraitsT>::run() {
  SPDLOG_INFO("LP NAME {} BASIS SIZE {}, COLUMN PIVOT RULE {}",
              _simplexTableau._initialProgram.getName(),
              _simplexTableau._rowInfos.size(),
              primalSimplexColumnPivotRuleToStr(_primalSimplexColumnPivotRule));
  auto artLpOptStats = primalSimplex().run(
      "PRIMAL_PHASE_ONE", PrintSimplexOptSummary::YES, PrimalPhase::ONE);

  if (_simplexTableau._objectiveValue >
      NumericalTraitsT::OBJECTIVE_MONOTONICITY_TOLERANCE) {
    SPDLOG_WARN(
        "PROGRAM WITH ARTIFICIAL VARIABLE HAS OPTIMUM {} GREATER THAN 0 - "
        "INITIAL PROGRAM IS INFEASIBLE",
        _simplexTableau._objectiveValue);
    artLpOptStats._phaseOneSucceeded = false;
    return artLpOptStats;
  }

  SimplexTableauResizer simplexTableauResizer(_simplexTableau,
                                              _reinversionManager);
  artLpOptStats._phaseOneSucceeded =
      simplexTableauResizer.removeArtificialVariablesFromProgram();
  return artLpOptStats;
}

template <typename T, typename SimplexTraitsT>
void PrimalSimplexPhaseOne<T, SimplexTraitsT>::makeRightHandSidesNonNegative() {
  SPDLOG_INFO("Making RHS non-negative");
  const int firstRowIdx =
      ((int)_simplexTableau._variableInfos[0]._isObjectiveVar);
  const int firstVarIdx = firstRowIdx;
  for (int rowIdx = firstRowIdx; rowIdx < _simplexTableau._rowInfos.size();
       ++rowIdx) {
    typename SimplexTraitsT::CurrentAdder adder;
    for (int variableIdx = firstVarIdx;
         variableIdx < _simplexTableau._variableInfos.size(); ++variableIdx) {
      adder.addValue(
          _simplexTableau._variableLowerBounds[variableIdx].value_or(0.0) *
          _simplexTableau._constraintMatrix[rowIdx][variableIdx]);
    }

    const T diff = NumericalTraitsT::add(
        _simplexTableau._rightHandSides[rowIdx], -adder.currentSum());

    if (diff < 0.0) {
      for (int variableIdx = firstVarIdx;
           variableIdx < _simplexTableau._variableInfos.size(); ++variableIdx)
        _simplexTableau._constraintMatrix[rowIdx][variableIdx] =
            -_simplexTableau._constraintMatrix[rowIdx][variableIdx];
      _simplexTableau._rightHandSides[rowIdx] =
          -_simplexTableau._rightHandSides[rowIdx];
      _simplexTableau._initialRightHandSides[rowIdx] =
          -_simplexTableau._initialRightHandSides[rowIdx];
    }
  }
}

template class PrimalSimplexPhaseOne<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class PrimalSimplexPhaseOne<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class PrimalSimplexPhaseOne<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class PrimalSimplexPhaseOne<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
