#ifndef GMISOLVER_PRIMALSIMPLEXGOMORY_H
#define GMISOLVER_PRIMALSIMPLEXGOMORY_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/IPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

template <typename T> class LinearProgram;
template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class LexicographicOptimizer;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT> class PrimalSimplex;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class PrimalSimplexGomory {
public:
  PrimalSimplexGomory(
      const LinearProgram<T> &primalLinearProgram,
      SimplexTableau<T, SimplexTraitsT> &dualSimplexTableau,
      ReinversionManager<T, SimplexTraitsT> &reinversionManager,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const ValidateSimplexOption validateSimplexOption,
      const SlackCutRemovalCondition slackCutRemovalCondition,
      const LexicographicReoptType lexicographicReoptType,
      const int cutRoundLimit);

  std::string type() const;

  IPOptStatistics<T> run(const LPOptimizationType lpOptimizationType,
                         const GomoryCutChoosingRule gomoryCutChoosingRule);

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  LPRelaxationStatistics<T> runImpl(const int relaxationNo);

  PrimalSimplex<T, SimplexTraitsT> primalSimplex() const;
  LexicographicOptimizer<T, SimplexTraitsT> lexicographicOptimizer() const;

  std::vector<int> collectFractionalDualCoordinates(
      const GomoryCutChoosingRule gomoryCutChoosingRule) const;

  void addCutColumns(const int relaxationNo,
                     const std::vector<int> &fractionalDualCoordinates) const;

  const LinearProgram<T> &_primalLinearProgram;
  SimplexTableau<T, SimplexTraitsT> &_dualSimplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
  const SlackCutRemovalCondition _slackCutRemovalCondition;
  const LexicographicReoptType _lexicographicReoptType;
  const int _cutRoundLimit;
};

#endif // GMISOLVER_PRIMALSIMPLEXGOMORY_H
