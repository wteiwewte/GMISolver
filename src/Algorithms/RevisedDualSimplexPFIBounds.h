#ifndef GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H
#define GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class RevisedDualSimplexPFIBounds {
public:
  RevisedDualSimplexPFIBounds(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      const DualSimplexRowPivotRule dualSimplexRowPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency);

  std::string name() const;

  void run();
  bool runOneIteration();

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  std::optional<int> chooseRow();
  std::optional<int> chooseRowFirstEligible();
  std::optional<int> chooseRowBiggestViolation();
  std::optional<int>
  chooseEnteringColumnIdx(const int pivotRowIdx, const VectorT &pivotRow,
                          const bool isPivotRowUnderLowerBound);

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
};

#endif // GMISOLVER_REVISEDDUALSIMPLEXPFIBOUNDS_H
