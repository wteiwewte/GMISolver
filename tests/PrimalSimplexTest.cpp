#include "Algorithms/ReinversionManager.h"
#include "Algorithms/RevisedDualSimplexPFIBounds.h"
#include "Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"
#include "tests/CommonDefs.h"
#include "tests/LPTestBase.h"

#include <filesystem>

#include <gtest/gtest.h>

template <typename T> struct PrimalSimplexOutput {
  LPOptStatistics<T> _phaseOneLpOptStats;
  std::optional<LPOptStatistics<T>> _phaseTwoLpOptStats;
};

template <typename T, typename SimplexTraitsT>
PrimalSimplexOutput<T>
runPrimalSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::PRIMAL,
      absl::GetFlag(FLAGS_use_product_form_of_inverse));
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>
      revisedPrimalSimplexPfiBounds(
          simplexTableau, reinversionManager,
          PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
          absl::GetFlag(FLAGS_obj_value_logging_frequency),
          absl::GetFlag(FLAGS_validate_simplex_option));
  auto phaseOneLpOptStats = revisedPrimalSimplexPfiBounds.runPhaseOne();
  if (!phaseOneLpOptStats._phaseOneSucceeded) {
    SPDLOG_WARN("PHASE ONE OF {} ALGORITHM FAILED",
                revisedPrimalSimplexPfiBounds.type());
    return {._phaseOneLpOptStats = phaseOneLpOptStats,
            ._phaseTwoLpOptStats = std::nullopt};
  }
  return {._phaseOneLpOptStats = phaseOneLpOptStats,
          ._phaseTwoLpOptStats = revisedPrimalSimplexPfiBounds.runPhaseTwo()};
}

template <typename T>
class PrimalSimplexTest : public LPTestBase<T>, public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_reinversion_frequency, 60);
  }
};

TYPED_TEST_SUITE_P(PrimalSimplexTest);
TYPED_TEST_P(PrimalSimplexTest, runPrimalSimplexAndCompareWithGurobi) {
  using FloatingPointT = std::tuple_element_t<0, typename TypeParam::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeParam::types>;
  constexpr auto PRIMAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/primal_simplex_working_instances";
  //      "../../tests/dual_simplex_working_instances";

  constexpr size_t PRIMAL_SIMPLEX_BASIS_SIZE_LIMIT = 200;
  const LPOptimizationType lpOptimizationType{
      LPOptimizationType::LINEAR_RELAXATION};
  this->solveAndCompareInstancesFromSets(
      PRIMAL_SIMPLEX_TEST_DIR_PATH, PRIMAL_SIMPLEX_BASIS_SIZE_LIMIT,
      lpOptimizationType,
      [&](const auto &linearProgram,
          const std::filesystem::path &modelFileMpsPath) {
        const auto primalSimplexOutput =
            runPrimalSimplexWithImplicitBounds<FloatingPointT, SimplexTraitsT>(
                linearProgram);
        const auto gurobiLPOptStats =
            GurobiOptimizer("", modelFileMpsPath)
                .optimize<FloatingPointT>(lpOptimizationType);
        if (primalSimplexOutput._phaseOneLpOptStats._optResult ==
            LPOptimizationResult::INFEASIBLE) {
          const std::set<LPOptimizationResult> infeasibleResults{
              LPOptimizationResult::INFEASIBLE,
              LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
          ASSERT_TRUE(infeasibleResults.contains(gurobiLPOptStats._optResult));
        } else {
          ASSERT_TRUE(primalSimplexOutput._phaseTwoLpOptStats.has_value());
          this->compareWithGurobi(*primalSimplexOutput._phaseTwoLpOptStats,
                                  gurobiLPOptStats);
        }
      },
      false);
}

REGISTER_TYPED_TEST_SUITE_P(PrimalSimplexTest,
                            runPrimalSimplexAndCompareWithGurobi);

using PrimalSimplexTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
// using PrimalSimplexTypes = ::testing::Types<
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::NORMAL>>,
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::SPARSE>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::SPARSE>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(PrimalSimplexTestSuite, PrimalSimplexTest,
                               PrimalSimplexTypes);
