#ifndef GMISOLVER_PHASEONEUTILITIES_H
#define GMISOLVER_PHASEONEUTILITIES_H

#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/SimplexBasisData.h"
#include "src/Util/LexReoptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT> class PrimalSimplex;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class PhaseOneUtilities {
public:
  PhaseOneUtilities(SimplexTableau<T, SimplexTraitsT> &simplexTableau);

  void addArtificialVariables(const SimplexType simplexType);
  std::vector<T> artificialObjective() const;
  std::optional<SimplexBasisData>
  createBasisFromArtificialVars(const SimplexType simplexType) const;

private:
  bool shouldMapZerothVarToZerothRowImpl(const SimplexType simplexType) const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
};

#endif // GMISOLVER_PHASEONEUTILITIES_H
