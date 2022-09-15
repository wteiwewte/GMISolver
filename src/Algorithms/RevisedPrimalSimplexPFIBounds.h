#ifndef GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
#define GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/SimplexTraits.h"

#include <optional>

template <typename T, typename SimplexTraitsT> class SimplexTableau;

template <typename T,
          typename SimplexTraitsT = SimplexTraits<T>>
class RevisedPrimalSimplexPFIBounds {
public:
  RevisedPrimalSimplexPFIBounds(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const int32_t objValueLoggingFrequency,
      const int32_t reinversionFrequency);

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
  void changeTableau(const PivotRowData<T> &pivotRowData,
                     const int enteringColumnIdx,
                     const std::vector<T> &enteringColumn);
  void moveVarFromOneBoundToAnother(const PivotRowData<T> &pivotRowData,
                                    const int enteringColumnIdx,
                                    const std::vector<T> &enteringColumn);

  bool removeArtificialVariablesFromBasis();
  void removeArtificialVariablesFromProgram();
  void removeRows(const std::vector<bool> &shouldRowBeRemoved);

  void setInitialObjective();
  void calculateDual();

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const int32_t _objValueLoggingFrequency;
  const int32_t _reinversionFrequency;
};

#endif // GMISOLVER_REVISEDPRIMALSIMPLEXPFIBOUNDS_H
