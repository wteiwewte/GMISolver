#ifndef GMISOLVER_DUALSIMPLEXGOMORYWITHPRIMALCUTS_H
#define GMISOLVER_DUALSIMPLEXGOMORYWITHPRIMALCUTS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class RevisedDualSimplexPFIBounds;
template <typename T, typename SimplexTraitsT> class RevisedPrimalSimplexPFIBounds;


template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class DualSimplexGomoryWithPrimalCuts {
public:
  DualSimplexGomoryWithPrimalCuts(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const DualSimplexRowPivotRule dualSimplexRowPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency);

  std::string type() const;

  void run(LPOptStatisticsVec<T>& lpOptStatisticsVec);
private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  RevisedDualSimplexPFIBounds<T, SimplexTraitsT> dualSimplex() const;
  RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT> primalSimplex() const;

  void addCutsForFractionalVars() const;
//  bool

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
};

#endif // GMISOLVER_DUALSIMPLEXGOMORYWITHPRIMALCUTS_H
