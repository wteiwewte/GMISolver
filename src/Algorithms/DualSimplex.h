#ifndef GMISOLVER_DUALSIMPLEX_H
#define GMISOLVER_DUALSIMPLEX_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class DualSimplex {
public:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  DualSimplex(SimplexTableau<T, SimplexTraitsT> &simplexTableau,
              ReinversionManager<T, SimplexTraitsT> &reinversionManager,
              const DualSimplexRowPivotRule dualSimplexRowPivotRule,
              const int32_t objValueLoggingFrequency,
              const ValidateSimplexOption validateSimplexOption);

  std::string type() const;

  LPOptStatistics<T> run(const std::string &lpNameSuffix,
                         const PrintSimplexOptSummary printSimplexOptSummary,
                         const DualPhase dualPhase);
  bool runOneIteration();
  bool pivot(const int pivotRowIdx,
             const std::optional<int> customEnteringColumnIdx,
             const bool isPivotRowUnderLowerBound);

private:
  bool isPivotRowUnderLowerBound(const int pivotRowIdx) const;

  void tryLogObjValue(const int iterCount);
  bool tryReinversion(const int iterCount, const DualPhase dualPhase,
                      const LPOptStatistics<T> &lpOptStatistics);
  bool tryValidateIteration(const int iterCount,
                            const LPOptStatistics<T> &lpOptStatistics);
  void tryValidateOptimalSolutions(const LPOptStatistics<T> &lpOptStatistics,
                                   const DualPhase dualPhase);
  bool checkIterationLimit(const int iterCount);
  bool checkObjectiveProgress(const LPOptStatistics<T> &lpOptStatistics);

  std::optional<int> chooseRow();
  std::optional<int> chooseRowFirstEligible();
  std::optional<int> chooseRowBiggestViolation();
  std::optional<int>
  chooseEnteringColumnIdx(const int pivotRowIdx, const VectorT &pivotRow,
                          const bool isPivotRowUnderLowerBound);

  bool isColumnAllowedToEnterBasis(const int colIdx) const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
};

#endif // GMISOLVER_DUALSIMPLEX_H
