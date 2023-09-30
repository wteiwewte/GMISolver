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

  bool removeArtificialVariablesFromProgram();
  void removeRows(const std::vector<bool> &shouldRowBeRemoved);
  void removeVariables(const std::vector<bool> &shouldVarBeRemoved);

private:
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  std::vector<bool> moveArtificialVariablesOutOfBasis();

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
};

#endif // GMISOLVER_SIMPLEXTABLEAURESIZER_H
