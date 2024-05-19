#ifndef GMISOLVER_DUALSIMPLEXPHASEONE_H
#define GMISOLVER_DUALSIMPLEXPHASEONE_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT> class DualSimplex;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class DualSimplexPhaseOne {
public:
  DualSimplexPhaseOne(SimplexTableau<T, SimplexTraitsT> &simplexTableau,
                      ReinversionManager<T, SimplexTraitsT> &reinversionManager,
                      const DualSimplexRowPivotRule dualSimplexRowPivotRule,
                      const int32_t objValueLoggingFrequency,
                      const ValidateSimplexOption validateSimplexOption);

  std::string type() const;

  LPOptStatistics<T> run();

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  DualSimplex<T, SimplexTraitsT> dualSimplex() const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;

  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
};

#endif // GMISOLVER_DUALSIMPLEXPHASEONE_H
