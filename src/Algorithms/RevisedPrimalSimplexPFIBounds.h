#ifndef GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
#define GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class RevisedPrimalSimplexPFIBounds {
public:
  RevisedPrimalSimplexPFIBounds(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency);

  std::string type() const;

  LPOptStatistics<T> runPhaseOne();
  LPOptStatistics<T> runPhaseTwo();
  void lexicographicReoptimization(const bool minimize, const std::string& lexOptId, LPOptStatisticsVec<T>& lpOptStatisticsVec);

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  LPOptStatistics<T> runImpl(const std::string& lpNameSuffix);
  void tryLogObjValue(const int iterCount);
  bool tryReinversion(const int iterCount);
  bool checkIterationLimit(const int iterCount);
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

  std::vector<T> singleVarObjective(const int varIdx, const bool minimize);
  void fixNonBasicVariables(int &varsFixedCount);
  void unfixAllVariables();

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
};

#endif // GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
