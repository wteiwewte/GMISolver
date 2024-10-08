#ifndef GMISOLVER_PRIMALSIMPLEXCUTS_H
#define GMISOLVER_PRIMALSIMPLEXCUTS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/IPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

template <typename T> class LinearProgram;
template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class LexicographicOptimizer;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT> class PrimalSimplex;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class PrimalSimplexCuts {
public:
  PrimalSimplexCuts(
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
                     const std::vector<int> &fractionalDualIndices,
                     IPOptStatistics<T> &ipOptStatistics);

  std::vector<T> getIthColumnOfBasisInverse(const int dualIdx) const;
  std::vector<T> getRVec(const std::vector<T> &ithColumnOfBasisInverse) const;
  std::vector<T> getBWithTildeVec(const int dualIdx,
                                  const std::vector<T> &rVec) const;
  T computeYBWithTildeProduct(const int dualIdx,
                              const std::vector<T> &rVec) const;

  std::string newCutVarLabel(const int relaxationNo, const int dualIdx) const;
  void addNewVar(const int relaxationNo, const int dualIdx,
                 const T yBWithTildeProduct);

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

#endif // GMISOLVER_PRIMALSIMPLEXCUTS_H
