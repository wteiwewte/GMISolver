#include "src/Algorithms/PrimalSimplexGomory.h"

#include "src/Algorithms/LexicographicOptimizer.h"
#include "src/Algorithms/PrimalSimplex.h"
#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexTableauResizer.h"
#include "src/Util/Time.h"

template <typename T, typename SimplexTraitsT>
PrimalSimplexGomory<T, SimplexTraitsT>::PrimalSimplexGomory(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption,
    const SlackCutRemovalCondition slackCutRemovalCondition,
    const LexicographicReoptType lexicographicReoptType,
    const int cutRoundLimit)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption),
      _slackCutRemovalCondition(slackCutRemovalCondition),
      _lexicographicReoptType(lexicographicReoptType),
      _cutRoundLimit(cutRoundLimit) {}

template <typename T, typename SimplexTraitsT>
std::string PrimalSimplexGomory<T, SimplexTraitsT>::type() const {
  return fmt::format(
      "PRIMAL SIMPLEX GOMORY WITH DUAL CUTS ({}, {})",
      simplexTableauTypeToStr(_simplexTableau._simplexTableauType),
      matrixRepresentationTypeToStr(SimplexTraitsT::matrixRepresentationType));
}

template <typename T, typename SimplexTraitsT>
IPOptStatistics<T> PrimalSimplexGomory<T, SimplexTraitsT>::run(
    const LPOptimizationType lpOptimizationType,
    const GomoryCutChoosingRule gomoryCutChoosingRule) {
  SPDLOG_INFO("PRIMAL GOMORY WITH {} LEXICOGRAPHIC REOPTIMIZATION",
              lexicographicReoptTypeToStr(_lexicographicReoptType));
  IPOptStatistics<T> ipOptStatistics{
      ._lpName = _simplexTableau.getName(),
      ._algorithmType = type(),
      ._reinversionFrequency = _reinversionManager.reinversionFrequency()};
  _simplexTableau.setObjective(_simplexTableau._initialProgram.getObjective());

  ipOptStatistics._elapsedTimeSec = executeAndMeasureTime([&] {
    int relaxationNo = 1;
    ipOptStatistics._lpRelaxationStats.emplace_back() = runImpl(relaxationNo);

    if (lpOptimizationType == LPOptimizationType::LINEAR_RELAXATION)
      return;

    if (_simplexTableau._simplexTableauType != SimplexTableauType::FULL) {
      SPDLOG_ERROR(
          "GOMORY CUTS NOT SUPPORTED FOR SIMPLEX TABLEAU TYPE {}",
          simplexTableauTypeToStr(_simplexTableau._simplexTableauType));
      return;
    }
  });

  ipOptStatistics._optimalValue = ipOptStatistics._lpRelaxationStats.back()
                                      ._lexicographicReoptStats._optimalValue;
  ipOptStatistics._optimalSolution = _simplexTableau._x;
  ipOptStatistics._optResult = _simplexTableau._result;

  return ipOptStatistics;
}

template <typename T, typename SimplexTraitsT>
LPRelaxationStatistics<T>
PrimalSimplexGomory<T, SimplexTraitsT>::runImpl(const int relaxationNo) {
  const auto relaxationId = [&relaxationNo] {
    return fmt::format("{}TH_RELAX", relaxationNo);
  };

  LPRelaxationStatistics<T> relaxationStats;
  relaxationStats._relaxationOptStats = primalSimplex().run(
      relaxationId(), PrintSimplexOptSummary::YES, PrimalPhase::TWO);
  _simplexTableau.calculateSolution();
  _simplexTableau.calculateCurrentObjectiveValue();
  SPDLOG_DEBUG(_simplexTableau.toStringObjectiveValue());
  SPDLOG_DEBUG(_simplexTableau.toStringSolution());

  relaxationStats._lexicographicReoptStats =
      lexicographicOptimizer().run(relaxationId());
  SPDLOG_DEBUG(_simplexTableau.toStringObjectiveValue());
  SPDLOG_DEBUG(_simplexTableau.toStringSolution());
  return relaxationStats;
}

template <typename T, typename SimplexTraitsT>
PrimalSimplex<T, SimplexTraitsT>
PrimalSimplexGomory<T, SimplexTraitsT>::primalSimplex() const {
  return PrimalSimplex<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _primalSimplexColumnPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
LexicographicOptimizer<T, SimplexTraitsT>
PrimalSimplexGomory<T, SimplexTraitsT>::lexicographicOptimizer() const {
  return LexicographicOptimizer<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _primalSimplexColumnPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption,
      _lexicographicReoptType);
}

template class PrimalSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class PrimalSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;