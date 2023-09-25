#include "src/Algorithms/DualSimplexGomory.h"

#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"

template <typename T, typename SimplexTraitsT>
DualSimplexGomory<T, SimplexTraitsT>::DualSimplexGomory(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const DualSimplexRowPivotRule dualSimplexRowPivotRule,
    const int32_t objValueLoggingFrequency, const int32_t reinversionFrequency,
    const ValidateSimplex validateSimplex)
    : _simplexTableau(simplexTableau),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _dualSimplexRowPivotRule(dualSimplexRowPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency),
      _validateSimplex(validateSimplex) {}

template <typename T, typename SimplexTraitsT>
std::string DualSimplexGomory<T, SimplexTraitsT>::type() const {
  return "DUAL SIMPLEX GOMORY WITH PRIMAL CUTS (" +
         std::string(SimplexTraitsT::useSparseRepresentationValue ? "SPARSE"
                                                                  : "NORMAL") +
         ')';
}

template <typename T, typename SimplexTraitsT>
IPOptStatistics<T> DualSimplexGomory<T, SimplexTraitsT>::run(
    const LexicographicReoptType lexicographicReoptType) {
  IPOptStatistics<T> ipOptStatistics;
  int relaxationCount = 1;
  const auto relaxationId = [&relaxationCount] {
    return fmt::format("{}TH_RELAX", relaxationCount);
  };

  auto &currentRelaxationStats =
      ipOptStatistics._lpRelaxationStats.emplace_back();
  currentRelaxationStats._relaxationOptStats =
      dualSimplex().run(relaxationId());
  SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());

  currentRelaxationStats._lexicographicReoptStats =
      primalSimplex().lexicographicReoptimization(lexicographicReoptType,
                                                  relaxationId());
  SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());

  //  while (true)
  //  {
  //
  //  }
  return ipOptStatistics;
}

template <typename T, typename SimplexTraitsT>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>
DualSimplexGomory<T, SimplexTraitsT>::dualSimplex() const {
  return RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
      _simplexTableau, _dualSimplexRowPivotRule, _objValueLoggingFrequency,
      _reinversionFrequency, _validateSimplex);
}

template <typename T, typename SimplexTraitsT>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>
DualSimplexGomory<T, SimplexTraitsT>::primalSimplex() const {
  return RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>(
      _simplexTableau, _primalSimplexColumnPivotRule, _objValueLoggingFrequency,
      _reinversionFrequency, _validateSimplex);
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomory<T, SimplexTraitsT>::addCutsForFractionalVars() const {
  return;
}

template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;