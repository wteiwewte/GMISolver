#ifndef GMISOLVER_LEXICOGRAPHICOPTIMIZER_H
#define GMISOLVER_LEXICOGRAPHICOPTIMIZER_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT>
class RevisedPrimalSimplexPFIBounds;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class LexicographicOptimizer {
public:
  LexicographicOptimizer(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency,
      const ValidateSimplex validateSimplex);

  std::string type() const;

  LexReoptStatistics<T> run(const LexicographicReoptType lexicographicReoptType,
                            const std::string &lexOptId,
                            const bool saveSolution = false);

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT> primalSimplex() const;

  std::vector<T>
  singleVarObjective(const int varIdx,
                     const LexicographicReoptType lexicographicReoptType);
  void fixNonBasicVariables(int &varsFixedCount);
  void unfixAllVariables();

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
  const ValidateSimplex _validateSimplex;
};

#endif // GMISOLVER_LEXICOGRAPHICOPTIMIZER_H
