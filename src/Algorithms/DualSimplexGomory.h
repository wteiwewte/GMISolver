#ifndef GMISOLVER_DUALSIMPLEXGOMORY_H
#define GMISOLVER_DUALSIMPLEXGOMORY_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/IPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class LexicographicOptimizer;
template <typename T, typename SimplexTraitsT>
class RevisedDualSimplexPFIBounds;
template <typename T, typename SimplexTraitsT>
class RevisedPrimalSimplexPFIBounds;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class DualSimplexGomory {
public:
  DualSimplexGomory(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const DualSimplexRowPivotRule dualSimplexRowPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency,
      const ValidateSimplexOption validateSimplexOption);

  std::string type() const;

  IPOptStatistics<T> run(const LexicographicReoptType lexicographicReoptType,
                         const LPOptimizationType lpOptimizationType,
                         const GomoryCutChoosingRule gomoryCutChoosingRule);

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  LPRelaxationStatistics<T>
  runImpl(const int relaxationNo,
          const LexicographicReoptType lexicographicReoptType);

  RevisedDualSimplexPFIBounds<T, SimplexTraitsT> dualSimplex() const;
  LexicographicOptimizer<T, SimplexTraitsT> lexicographicOptimizer() const;

  void checkIfNonBasicVarsAreIntegral() const;
  std::vector<int> collectFractionalBasisRowIndices(
      const GomoryCutChoosingRule gomoryCutChoosingRule) const;

  void addCutRows(const int relaxationNo,
                  const std::vector<int> &fractionalBasisVarsRowIndices) const;
  void
  addSlackVars(const int relaxationNo,
               const std::vector<int> &fractionalBasisVarsRowIndices) const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
  const ValidateSimplexOption _validateSimplexOption;
};

#endif // GMISOLVER_DUALSIMPLEXGOMORY_H
