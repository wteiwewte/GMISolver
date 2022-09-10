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
  RevisedPrimalSimplexPFIBounds(SimplexTableau<T> &simplexTableau,
                                         const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule);

  void run();
  bool runPhaseOne();
  void runPhaseTwo();
  void runImpl();
  bool runOneIteration();

private:
  std::optional<int> chooseEnteringColumn();
  std::optional<int> chooseEnteringColumnFirstEligible();
  std::optional<int> chooseEnteringColumnBiggestAbsReducedCost();

  std::optional<PivotRowData<T>>
  chooseRowIdx(const int enteringColumnIdx,
               const std::vector<T> &enteringColumn);
  void changeTableau(const PivotRowData<T> &pivotRowData, const int enteringColumnIdx,
             const std::vector<T> &enteringColumn);
  void moveVarFromOneBoundToAnother(const PivotRowData<T> &pivotRowData, const int enteringColumnIdx,
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
  PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
};

#endif // GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
