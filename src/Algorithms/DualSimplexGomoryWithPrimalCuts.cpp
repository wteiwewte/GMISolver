#include "src/Algorithms/DualSimplexGomoryWithPrimalCuts.h"

#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"

template <typename T, typename SimplexTraitsT>
DualSimplexGomoryWithPrimalCuts<T, SimplexTraitsT>::
    DualSimplexGomoryWithPrimalCuts(
        SimplexTableau<T, SimplexTraitsT> &simplexTableau,
        const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
        const DualSimplexRowPivotRule dualSimplexRowPivotRule,
        const int32_t objValueLoggingFrequency,
        const int32_t reinversionFrequency, const bool validateSimplex)
    : _simplexTableau(simplexTableau),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _dualSimplexRowPivotRule(dualSimplexRowPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency),
      _validateSimplex(validateSimplex) {}

template <typename T, typename SimplexTraitsT>
std::string DualSimplexGomoryWithPrimalCuts<T, SimplexTraitsT>::type() const {
  return "DUAL SIMPLEX GOMORY WITH PRIMAL CUTS (" +
         std::string(SimplexTraitsT::useSparseRepresentationValue ? "SPARSE"
                                                                  : "NORMAL") +
         ')';
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomoryWithPrimalCuts<T, SimplexTraitsT>::run(
    LPOptStatisticsVec<T> &lpOptStatisticsVec) {
  int relaxationCount = 1;
  const auto relaxationId = [&relaxationCount] {
    return fmt::format("{}TH_RELAX", relaxationCount);
  };

  auto lpStatisticsFromDualSimplex = dualSimplex().run(relaxationId());
  lpOptStatisticsVec.push_back(std::move(lpStatisticsFromDualSimplex));
  SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());

  primalSimplex().lexicographicReoptimization(false, relaxationId(),
                                              lpOptStatisticsVec);
  SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());

  //  while (true)
  //  {
  //
  //  }
}

template <typename T, typename SimplexTraitsT>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>
DualSimplexGomoryWithPrimalCuts<T, SimplexTraitsT>::dualSimplex() const {
  return RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
      _simplexTableau, _dualSimplexRowPivotRule, _objValueLoggingFrequency,
      _reinversionFrequency, _validateSimplex);
}

template <typename T, typename SimplexTraitsT>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>
DualSimplexGomoryWithPrimalCuts<T, SimplexTraitsT>::primalSimplex() const {
  return RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>(
      _simplexTableau, _primalSimplexColumnPivotRule, _objValueLoggingFrequency,
      _reinversionFrequency, _validateSimplex);
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomoryWithPrimalCuts<
    T, SimplexTraitsT>::addCutsForFractionalVars() const {
  return;
}

template class DualSimplexGomoryWithPrimalCuts<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexGomoryWithPrimalCuts<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;