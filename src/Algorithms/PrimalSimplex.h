#ifndef GMISOLVER_PRIMALSIMPLEX_H
#define GMISOLVER_PRIMALSIMPLEX_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/ComparisonTraits.h"

#include <optional>

template <typename T>
class SimplexTableau;

template <typename T, typename ComparisonTraitsT = ApproximateComparisonTraits<T>> class PrimalSimplex {
public:
  explicit PrimalSimplex(SimplexTableau<T>& simplexTableau);

  void runPhaseOne();

  void run();
  bool runOneIteration();
private:
  std::optional<int> chooseRowIdx(const int enteringColumnIdx);
  void pivot(const int rowIdx, const int enteringColumnIdx);

  void updateReducedCosts(const PivotData<T>& pivotData);
  void updateConstraintMatrixWithRHS(const PivotData<T>& pivotData);

  void removeArtificialVariablesFromBasis();
  void removeArtificialVariablesFromProgram();
  void removeRow(const int rowIdx);
  void setInitialObjective();

  SimplexTableau<T>& _simplexTableau;
};

#endif // GMISOLVER_PRIMALSIMPLEX_H
