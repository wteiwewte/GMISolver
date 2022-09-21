#ifndef GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDSSPARSE_H
#define GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDSSPARSE_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;

template <typename T,
          typename SimplexTraitsT = SimplexTraits<T>>
class RevisedPrimalSimplexPFIBoundsSparse {
public:
  RevisedPrimalSimplexPFIBoundsSparse(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency);

  void run();
  bool runPhaseOne();
  void runPhaseTwo();
  void runImpl();
  bool runOneIteration();

  void lexicographicReoptimization(const bool minimize);

private:
  std::optional<int> chooseEnteringColumn();
  std::optional<int> chooseEnteringColumnFirstEligible();
  std::optional<int> chooseEnteringColumnBiggestAbsReducedCost();

  std::optional<PivotRowData<T>>
  chooseRowIdx(const int enteringColumnIdx,
               const std::vector<T> &enteringColumn);
  void changeTableau(const PivotRowData<T> &pivotRowData,
                     const int enteringColumnIdx,
                     const SparseVector<T> &enteringColumn);
  void moveVarFromOneBoundToAnother(const PivotRowData<T> &pivotRowData,
                                    const int enteringColumnIdx,
                                    const std::vector<T> &enteringColumn);

  bool removeArtificialVariablesFromBasis();
  void removeArtificialVariablesFromProgram();
  void removeRows(const std::vector<bool> &shouldRowBeRemoved);

  std::vector<T> singleVarObjective(const int varIdx, const bool minimize);
  void fixNonBasicVariables(int& varsFixedCount);
  void unfixAllVariables();

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
};

#endif // GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDSSPARSE_H
