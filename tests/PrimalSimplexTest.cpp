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
PrimalSimplexOutput<T> runPrimalSimplexWithImplicitBounds(
    const LinearProgram<T> &linearProgram,
    const SimplexTableauType simplexTableauType) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::PRIMAL, simplexTableauType);
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
    absl::SetFlag(&FLAGS_simplex_tableau_types,
                  {SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE,
                   SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE,
                   SimplexTableauType::FULL});
    absl::SetFlag(&FLAGS_reinversion_frequency, 60);
  }

  void testCase(const std::string &primalSimplexTestDirPath,
                const size_t basisSizeLimit) {
    using FloatingPointT = std::tuple_element_t<0, typename T::types>;
    using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
    const LPOptimizationType lpOptimizationType{
        LPOptimizationType::LINEAR_RELAXATION};
    this->solveAndCompareInstancesFromSets(
        primalSimplexTestDirPath, basisSizeLimit, lpOptimizationType,
        [&](const auto &linearProgram,
            const SimplexTableauType simplexTableauType,
            const std::filesystem::path &modelFileMpsPath) {
          const auto primalSimplexOutput =
              runPrimalSimplexWithImplicitBounds<FloatingPointT,
                                                 SimplexTraitsT>(
                  linearProgram, simplexTableauType);
          const auto gurobiLPOptStats =
              GurobiOptimizer("", modelFileMpsPath)
                  .optimize<FloatingPointT>(lpOptimizationType);
          if (primalSimplexOutput._phaseOneLpOptStats._optResult ==
              LPOptimizationResult::INFEASIBLE) {
            const std::set<LPOptimizationResult> infeasibleResults{
                LPOptimizationResult::INFEASIBLE,
                LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
            ASSERT_TRUE(
                infeasibleResults.contains(gurobiLPOptStats._optResult));
          } else {
            ASSERT_TRUE(primalSimplexOutput._phaseTwoLpOptStats.has_value());
            this->compareWithGurobi(*primalSimplexOutput._phaseTwoLpOptStats,
                                    gurobiLPOptStats);
          }
        },
        false);
  }
};

TYPED_TEST_SUITE_P(PrimalSimplexTest);
TYPED_TEST_P(PrimalSimplexTest,
             runPrimalSimplexAndCompareWithGurobiBaseInstanceSet) {
  EXPECT_NO_FATAL_FAILURE(
      this->testCase("../../tests/primal_simplex_working_instances", 500));
}

TYPED_TEST_P(
    PrimalSimplexTest,
    runPrimalSimplexAndCompareWithGurobiNetlibInstanceSetWithFreeVars) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_simplex_tableau_types, {SimplexTableauType::FULL});
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/primal_simplex_netlib_instances_with_free_vars", 1000));
}

REGISTER_TYPED_TEST_SUITE_P(
    PrimalSimplexTest, runPrimalSimplexAndCompareWithGurobiBaseInstanceSet,
    runPrimalSimplexAndCompareWithGurobiNetlibInstanceSetWithFreeVars);

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
