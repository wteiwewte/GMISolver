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
  SPDLOG_INFO("Making RHS non-negative");
  _simplexTableau.makeRightHandSidesNonNegative();
  _simplexTableau.addArtificialVariables();
  _simplexTableau.initMatrixRepresentations();
  _simplexTableau.init(SimplexType::DUAL);

  _simplexTableau.calculateRHS();
  _simplexTableau.calculateSolution();
  _simplexTableau.calculateCurrentObjectiveValue();
  SPDLOG_TRACE("Simplex tableau with artificial variables");
  SPDLOG_TRACE(toString());
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
DualSimplex<T, SimplexTraitsT>
DualSimplexPhaseOne<T, SimplexTraitsT>::dualSimplex() const {
  return DualSimplex<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _dualSimplexRowPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T> DualSimplexPhaseOne<T, SimplexTraitsT>::run() {
  return LPOptStatistics<T>{._phaseOneSucceeded = true};
}

template class DualSimplexPhaseOne<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexPhaseOne<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class DualSimplexPhaseOne<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexPhaseOne<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
