#ifndef GMISOLVER_REVISEDPRIMALSIMPLEXPFI_H
#define GMISOLVER_REVISEDPRIMALSIMPLEXPFI_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/ComparisonTraits.h"

#include <optional>

template <typename T>
class SimplexTableau;

template <typename T, typename ComparisonTraitsT = ApproximateComparisonTraits<T>> class RevisedPrimalSimplexPFI {
public:
  explicit RevisedPrimalSimplexPFI(SimplexTableau<T>& simplexTableau);

  void runPhaseOne();

  void run();
  bool runOneIteration();
private:
  std::vector<T> computeEnteringColumn(const int enteringColumnIdx);
  std::vector<T> computePivotRow(const int rowIdx);
  std::optional<int> chooseRowIdx(const std::vector<T>& enteringColumn);
  void pivot(const int rowIdx, const int enteringColumnIdx, const std::vector<T>& enteringColumn, const std::vector<T>& pivotRow);

  void updateReducedCosts(const PivotData<T>& pivotData, const std::vector<T>& pivotRow);
  void updateInverseMatrixWithRHS(const PivotData<T>& pivotData, const std::vector<T>& enteringColumn);
  void calculateDual();



  void removeArtificialVariablesFromBasis();
  void removeArtificialVariablesFromProgram();
  void removeRow(const int rowIdx);
  void setInitialObjective();

  SimplexTableau<T>& _simplexTableau;
};

#endif // GMISOLVER_REVISEDPRIMALSIMPLEXPFI_H
