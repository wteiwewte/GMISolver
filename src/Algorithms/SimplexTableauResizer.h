#ifndef GMISOLVER_SIMPLEXTABLEAURESIZER_H
#define GMISOLVER_SIMPLEXTABLEAURESIZER_H

#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class SimplexTableauResizer {
public:
  SimplexTableauResizer(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      ReinversionManager<T, SimplexTraitsT> &reinversionManager);

  void removeArtificialVariablesFromBasis();
  void removeArtificialVariablesFromProgram();

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  void removeRows(const std::vector<bool> &shouldRowBeRemoved);

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
};

#endif // GMISOLVER_SIMPLEXTABLEAURESIZER_H
