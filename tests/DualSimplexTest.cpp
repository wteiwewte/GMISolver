#include "Algorithms/DualSimplex.h"
#include "Algorithms/DualSimplexPhaseOne.h"
#include "Algorithms/PrimalSimplex.h"
#include "Algorithms/ReinversionManager.h"
#include "Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"
#include "tests/CommonDefs.h"
#include "tests/LPTestBase.h"

#include <filesystem>

#include <gtest/gtest.h>

template <typename T, typename SimplexTraitsT>
SimplexOptimizationOutput<T>
runDualSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram,
                                 const SimplexTableauType simplexTableauType) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));

  DualSimplex<T, SimplexTraitsT> dualSimplex(
      simplexTableau, reinversionManager,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));

  DualSimplexPhaseOne<T, SimplexTraitsT> dualSimplexPhaseOne(
      simplexTableau, reinversionManager,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));

  auto phaseOneLpOptStats = dualSimplexPhaseOne.run();
  if (!phaseOneLpOptStats._phaseOneSucceeded) {
    SPDLOG_WARN("PHASE ONE OF {} ALGORITHM FAILED", dualSimplex.type());
    return {._phaseOneLpOptStats = phaseOneLpOptStats,
            ._phaseTwoLpOptStats = std::nullopt};
  }
  return {._phaseOneLpOptStats = phaseOneLpOptStats,
          ._phaseTwoLpOptStats =
              dualSimplex.run("", PrintSimplexOptSummary::YES, DualPhase::TWO)};
}

template <typename T>
class DualSimplexTest : public LPTestBase<T>, public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_simplex_tableau_types, {SimplexTableauType::FULL});
    absl::SetFlag(&FLAGS_extended_statistics, true);
  }

  void testCase(const std::string &dualSimplexTestDirPath,
                const size_t basisSizeLimit) {
    using FloatingPointT = std::tuple_element_t<0, typename T::types>;
    using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
    using UnderlyingOptStatsT = SimplexOptimizationOutput<FloatingPointT>;
    const LPOptimizationType lpOptimizationType{
        LPOptimizationType::LINEAR_RELAXATION};
    this->template solveAndCompareInstancesFromSets<UnderlyingOptStatsT>(
        dualSimplexTestDirPath, basisSizeLimit, lpOptimizationType,
        [&](const auto &linearProgram,
            const SimplexTableauType simplexTableauType,
            const std::filesystem::path &modelFileMpsPath,
            std::vector<OptimizationStats<UnderlyingOptStatsT>> &optStatsVec) {
          const auto dualSimplexLpOptStats =
              runDualSimplexWithImplicitBounds<FloatingPointT, SimplexTraitsT>(
                  linearProgram, simplexTableauType);
          const auto gurobiLPOptStats =
              GurobiOptimizer("", modelFileMpsPath)
                  .optimize<FloatingPointT>(
                      LPOptimizationType::LINEAR_RELAXATION);
          optStatsVec.push_back({{dualSimplexLpOptStats._phaseOneLpOptStats,
                                  dualSimplexLpOptStats._phaseTwoLpOptStats},
                                 gurobiLPOptStats});
          const bool isPrimalProgramDualFeasible =
              dualSimplexLpOptStats._phaseOneLpOptStats._phaseOneSucceeded &&
              dualSimplexLpOptStats._phaseTwoLpOptStats->_optResult ==
                  LPOptimizationResult::BOUNDED_AND_FEASIBLE;
          if (isPrimalProgramDualFeasible) {
            this->compareWithGurobi(*dualSimplexLpOptStats._phaseTwoLpOptStats,
                                    gurobiLPOptStats);
          } else {
            const std::set<LPOptimizationResult> gurobiInfeasibleResults{
                LPOptimizationResult::INFEASIBLE,
                LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
            EXPECT_TRUE(
                gurobiInfeasibleResults.contains(gurobiLPOptStats._optResult));
          }
        },
        false);
  }
};

TYPED_TEST_SUITE_P(DualSimplexTest);

TYPED_TEST_P(DualSimplexTest, runDualSimplexAndCompareWithGurobi) {
  EXPECT_NO_FATAL_FAILURE(
      this->testCase("../../tests/dual_simplex_working_instances", 1000));
}

TYPED_TEST_P(DualSimplexTest,
             runDualSimplexAndCompareWithGurobiOnInfeasibleNetlibInstances) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_simplex_tableau_types, {SimplexTableauType::FULL});
  EXPECT_NO_FATAL_FAILURE(
      this->testCase("../../tests/netlib_infeasible_instances", 1000));
}

REGISTER_TYPED_TEST_SUITE_P(
    DualSimplexTest, runDualSimplexAndCompareWithGurobi,
    runDualSimplexAndCompareWithGurobiOnInfeasibleNetlibInstances);

using DualSimplexTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>,
    TypeTuple<long double,
              SimplexTraits<long double, MatrixRepresentationType::NORMAL>>>;
//    TypeTuple<double, SimplexTraits<double,
//    MatrixRepresentationType::SPARSE>>, TypeTuple<long double,
//              SimplexTraits<long double, MatrixRepresentationType::SPARSE>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(DualSimplexTestSuite, DualSimplexTest,
                               DualSimplexTypes);
