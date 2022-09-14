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
  RevisedDualSimplexPFIBounds(
      SimplexTableau<T> &simplexTableau,
      const DualSimplexRowPivotRule dualSimplexRowPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency);

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
  void pivot(const std::vector<T> &pivotRow, const int pivotRowIdx,
             const std::vector<T> &enteringColumn, const int enteringColumnIdx,
             const bool isPivotRowUnderLowerBound);

  void calculateCurrentObjectiveValue();
  void calculateSolution();
  void initRHS();

  SimplexTableau<T> &_simplexTableau;
  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
};

#endif // GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H
