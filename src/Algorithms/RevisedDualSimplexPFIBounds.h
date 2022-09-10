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
  RevisedDualSimplexPFIBounds(SimplexTableau<T> &simplexTableau,
                                       const DualSimplexRowPivotRule dualSimplexRowPivotRule);

  void run();
  bool runOneIteration();

private:
  std::optional<int> chooseRow();
  std::optional<int> chooseRowFirstEligible();
  std::optional<int> chooseRowBiggestViolation();
  std::vector<T> computePivotRow(const int rowIdx);
  std::optional<int>
  chooseEnteringColumnIdx(const int pivotRowIdx, const std::vector<T> &pivotRow,
                          const bool isPivotRowUnderLowerBound);
  //  std::optional<int> chooseEnteringColumnIdxBiggestReducedCost();
  //  std::optional<PivotRowData<T>>
  //  chooseRowIdx(const int enteringColumnIdx,
  //               const std::vector<T> &enteringColumn);
  void pivot(const std::vector<T> &pivotRow, const int pivotRowIdx,
             const std::vector<T> &enteringColumn, const int enteringColumnIdx,
             const bool isPivotRowUnderLowerBound);
  //  void calculateDual();
  //  void reinversion();
  //
  //
  void calculateCurrentObjectiveValue();
  void calculateSolution();
  void initRHS();
  //  void calculateRHS();

  SimplexTableau<T> &_simplexTableau;
  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
};

#endif // GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H
