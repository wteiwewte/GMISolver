#ifndef GMISOLVER_PRIMALSIMPLEX_H
#define GMISOLVER_PRIMALSIMPLEX_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class PrimalSimplex {
public:
  PrimalSimplex(SimplexTableau<T, SimplexTraitsT> &simplexTableau,
                ReinversionManager<T, SimplexTraitsT> &reinversionManager,
                const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
                const int32_t objValueLoggingFrequency,
                const ValidateSimplexOption validateSimplexOption);

  std::string type() const;

  LPOptStatistics<T> runPhaseTwo();
  LPOptStatistics<T> run(const std::string &lpNameSuffix,
                         const PrintSimplexOptSummary printSimplexOptSummary =
                             PrintSimplexOptSummary::YES,
                         const PrimalPhase primalPhase = PrimalPhase::TWO);

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  void tryLogObjValue(const int iterCount);
  bool tryReinversion(const int iterCount,
                      const LPOptStatistics<T> &lpOptStatistics,
                      const PrimalPhase primalPhase);
  bool checkIterationLimit(const int iterCount);
  bool checkObjectiveProgress(const LPOptStatistics<T> &lpOptStatistics);
  bool tryValidateIteration(const int iterCount,
                            const LPOptStatistics<T> &lpOptStatistics,
                            const PrimalPhase primalPhase);
  void tryValidateOptimalSolutions(const LPOptStatistics<T> &lpOptStatistics,
                                   const PrimalPhase primalPhase);
  bool runOneIteration();

  std::optional<int> chooseEnteringColumn();
  std::optional<int> chooseEnteringColumnFirstEligible();
  std::optional<int> chooseEnteringColumnBiggestAbsReducedCost();

  std::optional<PivotRowData<T>> chooseRowIdx(const int enteringColumnIdx,
                                              const VectorT &enteringColumn);
  void changeTableau(const PivotRowData<T> &pivotRowData,
                     const int enteringColumnIdx,
                     const VectorT &enteringColumn);
  void moveVarFromOneBoundToAnother(const PivotRowData<T> &pivotRowData,
                                    const int enteringColumnIdx,
                                    const VectorT &enteringColumn);
  bool isColumnAllowedToEnterBasis(const int colIdx) const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;

  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
};

#endif // GMISOLVER_PRIMALSIMPLEX_H
