#ifndef GMISOLVER_DUALSIMPLEXGOMORY_H
#define GMISOLVER_DUALSIMPLEXGOMORY_H

#include "src/DataModel/CommonTypes.h"
#include "src/Util/IPOptStatistics.h"
#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT> class SimplexTableau;
template <typename T, typename SimplexTraitsT> class LexicographicOptimizer;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT>
class RevisedDualSimplexPFIBounds;
template <typename T, typename SimplexTraitsT>
class RevisedPrimalSimplexPFIBounds;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class DualSimplexGomory {
public:
  DualSimplexGomory(
      SimplexTableau<T, SimplexTraitsT> &simplexTableau,
      ReinversionManager<T, SimplexTraitsT> &reinversionManager,
      const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
      const DualSimplexRowPivotRule dualSimplexRowPivotRule,
      const int32_t objValueLoggingFrequency,
      const ValidateSimplexOption validateSimplexOption);

  std::string type() const;

  IPOptStatistics<T> run(const LexicographicReoptType lexicographicReoptType,
                         const LPOptimizationType lpOptimizationType,
                         const GomoryCutChoosingRule gomoryCutChoosingRule);

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  LPRelaxationStatistics<T>
  runImpl(const int relaxationNo,
          const LexicographicReoptType lexicographicReoptType);

  RevisedDualSimplexPFIBounds<T, SimplexTraitsT> dualSimplex() const;
  LexicographicOptimizer<T, SimplexTraitsT> lexicographicOptimizer() const;

  void checkIfNonBasicVarsAreIntegral() const;
  std::vector<int> collectFractionalBasisRowIndices(
      const GomoryCutChoosingRule gomoryCutChoosingRule) const;

  void addCutRows(const int relaxationNo,
                  const std::vector<int> &fractionalBasisVarsRowIndices) const;
  void
  addSlackVars(const int relaxationNo,
               const std::vector<int> &fractionalBasisVarsRowIndices) const;

  bool removeCutsInBasis() const;
  void removeCutFromBasis(const int basisRowIdxMappedToCutVar,
                          const int cutVarIdx, const int cutRowIdx) const;

  void removeCuts(const std::vector<bool> &shouldVarBeRemoved,
                  const std::vector<bool> &shouldRowBeRemoved) const;

  void
  debugLogOldAndNewBasis(const std::vector<int> &oldRowIdxToNewRowIdx,
                         const std::vector<int> &newRowToBasicColumnIdxMap,
                         const std::vector<bool> &shouldRowBeRemoved) const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  ReinversionManager<T, SimplexTraitsT> &_reinversionManager;
  const PrimalSimplexColumnPivotRule _primalSimplexColumnPivotRule;
  const DualSimplexRowPivotRule _dualSimplexRowPivotRule;
  const int32_t _objValueLoggingFrequency;
  const ValidateSimplexOption _validateSimplexOption;
  const bool _defineCutsInTermsOfOriginalVariables{true};
};

#endif // GMISOLVER_DUALSIMPLEXGOMORY_H
