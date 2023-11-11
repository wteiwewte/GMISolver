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

  void testCaseWithDual(const std::string &primalSimplexTestDirPath,
                        const size_t basisSizeLimit) {
    using FloatingPointT = std::tuple_element_t<0, typename T::types>;
    using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
    const LPOptimizationType lpOptimizationType{
        LPOptimizationType::LINEAR_RELAXATION};
    this->solveAndCompareInstancesFromSets(
        primalSimplexTestDirPath, basisSizeLimit, lpOptimizationType,
        [&](const auto &primalProgram,
            const SimplexTableauType simplexTableauType,
            const std::filesystem::path &modelFileMpsPath) {
          const auto primalProgramSimplexOutput =
              runPrimalSimplexWithImplicitBounds<FloatingPointT,
                                                 SimplexTraitsT>(
                  primalProgram, simplexTableauType);
          const auto dualProgram = primalProgram.dualProgram();
          ASSERT_TRUE(dualProgram.has_value());
          const auto dualProgramSimplexOutput =
              runPrimalSimplexWithImplicitBounds<FloatingPointT,
                                                 SimplexTraitsT>(
                  *dualProgram, simplexTableauType);
          const auto gurobiLPOptStats =
              GurobiOptimizer("", modelFileMpsPath)
                  .optimize<FloatingPointT>(lpOptimizationType);

          ASSERT_EQ(LPOptimizationResult::BOUNDED_AND_FEASIBLE,
                    primalProgramSimplexOutput._phaseOneLpOptStats._optResult);
          const bool isPrimalProgramInfeasible =
              primalProgramSimplexOutput._phaseOneLpOptStats._optimalValue >
              SimplexTraitsT::NumericalTraitsT::
                  OBJECTIVE_MONOTONICITY_TOLERANCE;
          if (isPrimalProgramInfeasible) {
            const std::set<LPOptimizationResult> gurobiInfeasibleResults{
                LPOptimizationResult::INFEASIBLE,
                LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
            EXPECT_TRUE(
                gurobiInfeasibleResults.contains(gurobiLPOptStats._optResult));
            EXPECT_EQ(LPOptimizationResult::BOUNDED_AND_FEASIBLE,
                      dualProgramSimplexOutput._phaseOneLpOptStats._optResult);
            ASSERT_TRUE(
                dualProgramSimplexOutput._phaseTwoLpOptStats.has_value());
            EXPECT_EQ(LPOptimizationResult::UNBOUNDED,
                      dualProgramSimplexOutput._phaseTwoLpOptStats->_optResult);
          } else {
            ASSERT_TRUE(
                primalProgramSimplexOutput._phaseTwoLpOptStats.has_value());
            switch (
                primalProgramSimplexOutput._phaseTwoLpOptStats->_optResult) {
            case LPOptimizationResult::UNBOUNDED: {
              const std::set<LPOptimizationResult> gurobiUnboundedResults{
                  LPOptimizationResult::UNBOUNDED,
                  LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
              EXPECT_TRUE(
                  gurobiUnboundedResults.contains(gurobiLPOptStats._optResult));
              EXPECT_EQ(
                  LPOptimizationResult::INFEASIBLE,
                  dualProgramSimplexOutput._phaseOneLpOptStats._optResult);
              ASSERT_FALSE(
                  dualProgramSimplexOutput._phaseTwoLpOptStats.has_value());
              break;
            }
            case LPOptimizationResult::BOUNDED_AND_FEASIBLE: {
              EXPECT_EQ(
                  LPOptimizationResult::BOUNDED_AND_FEASIBLE,
                  dualProgramSimplexOutput._phaseOneLpOptStats._optResult);
              ASSERT_TRUE(
                  dualProgramSimplexOutput._phaseTwoLpOptStats.has_value());
              EXPECT_EQ(
                  LPOptimizationResult::BOUNDED_AND_FEASIBLE,
                  dualProgramSimplexOutput._phaseTwoLpOptStats->_optResult);
              EXPECT_NO_FATAL_FAILURE(this->compareWithGurobiDual(
                  *primalProgramSimplexOutput._phaseTwoLpOptStats,
                  *dualProgramSimplexOutput._phaseTwoLpOptStats,
                  gurobiLPOptStats));
              break;
            }
            default: {
              SPDLOG_INFO("COULD NOT FINISH PHASE TWO SIMPLEX FOR PRIMAL "
                          "PROGRAM DUE TO {}",
                          lpOptimizationResultToStr(
                              primalProgramSimplexOutput._phaseTwoLpOptStats
                                  ->_optResult));
              EXPECT_TRUE(false);
              break;
            }
            }
          }
        },
        false);
  }

  void testCase(const std::string &primalSimplexTestDirPath,
                const size_t basisSizeLimit) {
    using FloatingPointT = std::tuple_element_t<0, typename T::types>;
    using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
    const LPOptimizationType lpOptimizationType{
        LPOptimizationType::LINEAR_RELAXATION};
    this->solveAndCompareInstancesFromSets(
        primalSimplexTestDirPath, basisSizeLimit, lpOptimizationType,
        [&](const auto &primalProgram,
            const SimplexTableauType simplexTableauType,
            const std::filesystem::path &modelFileMpsPath) {
          const auto primalProgramSimplexOutput =
              runPrimalSimplexWithImplicitBounds<FloatingPointT,
                                                 SimplexTraitsT>(
                  primalProgram, simplexTableauType);
          const auto gurobiLPOptStats =
              GurobiOptimizer("", modelFileMpsPath)
                  .optimize<FloatingPointT>(lpOptimizationType);
          ASSERT_EQ(LPOptimizationResult::BOUNDED_AND_FEASIBLE,
                    primalProgramSimplexOutput._phaseOneLpOptStats._optResult);
          const bool isPrimalProgramInfeasible =
              primalProgramSimplexOutput._phaseOneLpOptStats._optimalValue >
              SimplexTraitsT::NumericalTraitsT::
                  OBJECTIVE_MONOTONICITY_TOLERANCE;
          if (isPrimalProgramInfeasible) {
            const std::set<LPOptimizationResult> gurobiInfeasibleResults{
                LPOptimizationResult::INFEASIBLE,
                LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
            EXPECT_TRUE(
                gurobiInfeasibleResults.contains(gurobiLPOptStats._optResult));
          } else {
            ASSERT_TRUE(
                primalProgramSimplexOutput._phaseTwoLpOptStats.has_value());
            switch (
                primalProgramSimplexOutput._phaseTwoLpOptStats->_optResult) {
            case LPOptimizationResult::UNBOUNDED: {
              const std::set<LPOptimizationResult> gurobiUnboundedResults{
                  LPOptimizationResult::UNBOUNDED,
                  LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
              EXPECT_TRUE(
                  gurobiUnboundedResults.contains(gurobiLPOptStats._optResult));
              break;
            }
            case LPOptimizationResult::BOUNDED_AND_FEASIBLE: {
              EXPECT_NO_FATAL_FAILURE(this->compareWithGurobi(
                  *primalProgramSimplexOutput._phaseTwoLpOptStats,
                  gurobiLPOptStats));
              break;
            }
            default: {
              SPDLOG_INFO("COULD NOT FINISH PHASE TWO SIMPLEX FOR PRIMAL "
                          "PROGRAM DUE TO {}",
                          lpOptimizationResultToStr(
                              primalProgramSimplexOutput._phaseTwoLpOptStats
                                  ->_optResult));
              EXPECT_TRUE(false);
              break;
            }
            }
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
    runPrimalSimplexForPrimalAndDualAndCompareWithGurobiBaseInstanceSet) {
  EXPECT_NO_FATAL_FAILURE(this->testCaseWithDual(
      "../../tests/primal_simplex_working_instances", 200));
}

TYPED_TEST_P(
    PrimalSimplexTest,
    runPrimalSimplexForPrimalAndDualAndCompareWithGurobiBaseInstanceSetOnlyFull) {
  absl::SetFlag(&FLAGS_simplex_tableau_types, {SimplexTableauType::FULL});
  EXPECT_NO_FATAL_FAILURE(this->testCaseWithDual(
      "../../tests/primal_simplex_working_instances", 500));
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

TYPED_TEST_P(PrimalSimplexTest,
             runPrimalSimplexAndCompareWithGurobiNetlibInfeasibleInstanceSet) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_simplex_tableau_types, {SimplexTableauType::FULL});
  absl::SetFlag(&FLAGS_obj_value_logging_frequency, 4000);
  EXPECT_NO_FATAL_FAILURE(
      this->testCase("../../tests/netlib_infeasible_instances", 400));
}
TYPED_TEST_P(
    PrimalSimplexTest,
    runPrimalSimplexForPrimalAndDualAndCompareWithGurobiNetlibInfeasibleInstanceSet) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_simplex_tableau_types, {SimplexTableauType::FULL});
  absl::SetFlag(&FLAGS_obj_value_logging_frequency, 4000);
  EXPECT_NO_FATAL_FAILURE(
      this->testCaseWithDual("../../tests/netlib_infeasible_instances", 300));
}

REGISTER_TYPED_TEST_SUITE_P(
    PrimalSimplexTest, runPrimalSimplexAndCompareWithGurobiBaseInstanceSet,
    runPrimalSimplexForPrimalAndDualAndCompareWithGurobiBaseInstanceSet,
    runPrimalSimplexForPrimalAndDualAndCompareWithGurobiBaseInstanceSetOnlyFull,
    runPrimalSimplexAndCompareWithGurobiNetlibInstanceSetWithFreeVars,
    runPrimalSimplexAndCompareWithGurobiNetlibInfeasibleInstanceSet,
    runPrimalSimplexForPrimalAndDualAndCompareWithGurobiNetlibInfeasibleInstanceSet);

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
