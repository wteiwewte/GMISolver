#ifndef GMISOLVER_LEXICOGRAPHICOPTIMIZER_H
#define GMISOLVER_LEXICOGRAPHICOPTIMIZER_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LexReoptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <boost/dynamic_bitset.hpp>

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT> class PrimalSimplex;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class LexicographicOptimizer {
public:
  LexicographicOptimizer(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      ReinversionManager<T, SimplexTraitsT> &reinversionManager,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const ValidateSimplexOption validateSimplexOption,
      const LexicographicReoptType lexicographicReoptType);

  std::string type() const;

  LexReoptStatistics<T>
  run(const std::string &lexOptId,
      const SaveLexSolution saveLexSolution = SaveLexSolution::NO);

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  PrimalSimplex<T, SimplexTraitsT> primalSimplex() const;

  boost::dynamic_bitset<> getIsFixedBitset() const;
  std::vector<T> singleVarObjective(const int varIdx);
  void fixNonBasicVariables(int &varsFixedCount);
  void unfixVariables(const boost::dynamic_bitset<> &initialIsFixedBitset);

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
  const LexicographicReoptType _lexicographicReoptType;
};

#endif // GMISOLVER_LEXICOGRAPHICOPTIMIZER_H
