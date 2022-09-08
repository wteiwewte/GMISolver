#ifndef GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H
#define GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/ComparisonTraits.h"

#include <optional>

template <typename T> class SimplexTableau;

template <typename T,
          typename ComparisonTraitsT = ApproximateComparisonTraits<T>>
class RevisedDualSimplexPFIBounds {
public:
  explicit RevisedDualSimplexPFIBounds(SimplexTableau<T> &simplexTableau);

  void run();
  bool runOneIteration();

private:
  std::optional<int> chooseRowIdx();
  std::optional<int> chooseRowIdxBiggestViolation();
  std::vector<T> computePivotRow(const int rowIdx);
  std::optional<int>
  chooseEnteringColumnIdx(const int pivotRowIdx, const std::vector<T> &pivotRow,
                          const bool isPivotRowUnderLowerBound);
  //  std::optional<int> chooseEnteringColumnIdxBiggestReducedCost();
  std::vector<T> computeEnteringColumn(const int enteringColumnIdx);
  //  std::optional<PivotRowData<T>>
  //  chooseRowIdx(const int enteringColumnIdx,
  //               const std::vector<T> &enteringColumn);
  void pivot(const std::vector<T> &pivotRow, const int pivotRowIdx,
             const std::vector<T> &enteringColumn, const int enteringColumnIdx,
             const bool isPivotRowUnderLowerBound);
  void executePivot(const int rowIdx, const int enteringColumnIdx,
                    const std::vector<T> &enteringColumn,
                    const std::vector<T> &pivotRow);
  void updateReducedCosts(const PivotData<T> &pivotData,
                          const std::vector<T> &pivotRow);
  void updateInverseMatrixWithRHS(const PivotData<T> &pivotData,
                                  const std::vector<T> &enteringColumn);
  //  void calculateDual();
  //  void reinversion();
  //
  //
  void calculateCurrentObjectiveValue();
  void calculateSolution();
  void initRHS();
  //  void calculateRHS();

  SimplexTableau<T> &_simplexTableau;
};

#endif // GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H
