#ifndef GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
#define GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class ReinversionManager;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class RevisedPrimalSimplexPFIBounds {
public:
  RevisedPrimalSimplexPFIBounds(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      ReinversionManager<T, SimplexTraitsT> &reinversionManager,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const ValidateSimplexOption validateSimplexOption);

  std::string type() const;

  LPOptStatistics<T> runPhaseOne();
  LPOptStatistics<T> runPhaseTwo();
  LPOptStatistics<T> runImpl(const std::string &lpNameSuffix,
                             const bool printSummary = true);

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  void tryLogObjValue(const int iterCount);
  bool tryReinversion(const int iterCount,
                      const LPOptStatistics<T> &lpOptStatistics);
  bool checkIterationLimit(const int iterCount);
  bool tryValidateIteration(const int iterCount,
                            const LPOptStatistics<T> &lpOptStatistics);
  void tryValidateOptimalSolutions(const LPOptStatistics<T> &lpOptStatistics);
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

  void removeArtificialVariablesFromBasis();
  void removeArtificialVariablesFromProgram();
  void removeRows(const std::vector<bool> &shouldRowBeRemoved);

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;

  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
};

#endif // GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
