#ifndef GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
#define GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/ComparisonTraits.h"

#include <optional>

template <typename T> class SimplexTableau;

template <typename T,
          typename ComparisonTraitsT = ApproximateComparisonTraits<T>>
class RevisedPrimalSimplexPFIBounds {
public:
  explicit RevisedPrimalSimplexPFIBounds(SimplexTableau<T> &simplexTableau);

  void run();
  bool runPhaseOne();
  void runPhaseTwo();
  void runImpl();
  bool runOneIteration();

private:
  std::optional<int> chooseEnteringColumnIdx();
  std::optional<int> chooseEnteringColumnIdxBiggestReducedCost();
  std::vector<T> computeEnteringColumn(const int enteringColumnIdx);
  std::vector<T> computePivotRow(const int rowIdx);
  std::optional<PivotRowData<T>>
  chooseRowIdx(const int enteringColumnIdx,
               const std::vector<T> &enteringColumn);
  void pivot(const PivotRowData<T> &pivotRowData, const int enteringColumnIdx,
             const std::vector<T> &enteringColumn);
  void executePivot(const int rowIdx, const int enteringColumnIdx,
                    const std::vector<T> &enteringColumn,
                    const std::vector<T> &pivotRow);
  void updateReducedCosts(const PivotData<T> &pivotData,
                          const std::vector<T> &pivotRow);
  void updateInverseMatrixWithRHS(const PivotData<T> &pivotData,
                                  const std::vector<T> &enteringColumn);
  void calculateDual();
  void reinversion();
  //  void calculateRHS();
  //
  void removeArtificialVariablesFromBasis();
  void removeArtificialVariablesFromProgram();
  void removeRow(const int rowIdx);
  void setInitialObjective();
  void calculateCurrentObjectiveValue();
  void calculateSolution();
  void initRHS();
  void calculateRHS();
  //  void initSolution();

  SimplexTableau<T> &_simplexTableau;
};

#endif // GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
