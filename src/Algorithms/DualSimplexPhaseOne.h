#ifndef GMISOLVER_DUALSIMPLEXPHASEONE_H
#define GMISOLVER_DUALSIMPLEXPHASEONE_H

#include "src/Algorithms/PhaseOneUtilities.h"
#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

#include <boost/dynamic_bitset.hpp>

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

  void initBoundsForBoxedVars();
  void markBoxedVariablesAsNotEligible();
  void unmarkBoxedVariablesAsNotEligible();
  void changeBoundsRHSAndObj();
  void restoreBoundsRHSAndObj();
  void recalculateRHS();
  std::vector<T> auxiliaryObjective() const;

  DualSimplex<T, SimplexTraitsT> dualSimplex() const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
  PhaseOneUtilities<T, SimplexTraitsT> _phaseOneUtilities;

  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;

  std::vector<VariableInfo> _originalVariableInfos;
  std::vector<std::optional<T>> _originalVariableLowerBounds;
  std::vector<std::optional<T>> _originalVariableUpperBounds;
  std::vector<T> _originalInitialRightHandSides;
};

#endif // GMISOLVER_DUALSIMPLEXPHASEONE_H
