#ifndef GMISOLVER_PRIMALSIMPLEXPHASEONE_H
#define GMISOLVER_PRIMALSIMPLEXPHASEONE_H

#include "src/Algorithms/PhaseOneUtilities.h"
#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT> class PrimalSimplex;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class PrimalSimplexPhaseOne {
public:
  PrimalSimplexPhaseOne(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      ReinversionManager<T, SimplexTraitsT> &reinversionManager,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const ValidateSimplexOption validateSimplexOption);

  std::string type() const;

  LPOptStatistics<T> run();

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  PrimalSimplex<T, SimplexTraitsT> primalSimplex() const;

  void makeRightHandSidesNonNegative();

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
  PhaseOneUtilities<T, SimplexTraitsT> _phaseOneUtilities;

  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
};

#endif // GMISOLVER_PRIMALSIMPLEXPHASEONE_H
