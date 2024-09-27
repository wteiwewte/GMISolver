#include "src/Algorithms/PrimalSimplexGomory.h"
#include "src/Algorithms/PrimalSimplexPhaseOne.h"
#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"
#include "tests/CommonDefs.h"
#include "tests/LPTestBase.h"

#include <filesystem>

#include <gtest/gtest.h>

template <typename T, typename SimplexTraitsT>
IPOptStatistics<T>
runPrimalSimplexGomory(const LinearProgram<T> &primalProgram,
                       const LinearProgram<T> &dualProgram,
                       const SimplexTableauType simplexTableauType,
                       const LexicographicReoptType lexicographicReoptType,
                       const LPOptimizationType lpOptimizationType,
                       const GomoryCutChoosingRule gomoryCutChoosingRule) {
  SimplexTableau<T, SimplexTraitsT> dualSimplexTableau(
      dualProgram, SimplexType::PRIMAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      dualSimplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));

  PrimalSimplexPhaseOne<T, SimplexTraitsT> primalSimplexPhaseOne(
      dualSimplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));
  auto phaseOneLpOptStats = primalSimplexPhaseOne.run();
  if (!phaseOneLpOptStats._phaseOneSucceeded) {
    SPDLOG_WARN("PHASE ONE OF {} ALGORITHM FAILED",
                primalSimplexPhaseOne.type());
    return {};
  }

  PrimalSimplexGomory<T, SimplexTraitsT> primalSimplexGomory(
      primalProgram, dualSimplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option),
      absl::GetFlag(FLAGS_slack_cut_removal_condition), lexicographicReoptType,
      absl::GetFlag(FLAGS_cut_round_limit));
  return primalSimplexGomory.run(lpOptimizationType, gomoryCutChoosingRule);
}

template <typename T>
class PrimalSimplexGomoryTest : public LPTestBase<T>, public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_simplex_tableau_types,
                  {SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE});
    absl::SetFlag(&FLAGS_slack_cut_removal_condition,
                  SlackCutRemovalCondition::ONLY_WHEN_SLACK_VAR_IS_POSITIVE);
    absl::SetFlag(&FLAGS_cut_round_limit, 2000);
    absl::SetFlag(&FLAGS_reinversion_frequency, 0);
  }

  void
  testCase(const std::string &PrimalSimplexGomoryTestDirPath,
           const size_t basisSizeLimit,
           const LPOptimizationType lpOptimizationType,
           const std::vector<LexicographicReoptType> lexicographicReoptTypes) {
    using FloatingPointT = std::tuple_element_t<0, typename T::types>;
    using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
    using UnderlyingOptStatsT = IPOptStatistics<FloatingPointT>;
    this->template solveAndCompareInstancesFromSets<UnderlyingOptStatsT>(
        PrimalSimplexGomoryTestDirPath, basisSizeLimit, lpOptimizationType,
        [&](const auto &primalProgram,
            const SimplexTableauType simplexTableauType,
            const std::filesystem::path &modelFileMpsPath,
            std::vector<OptimizationStats<UnderlyingOptStatsT>> &optStatsVec) {
          const auto dualProgram =
              primalProgram.dualProgram(AddObjectiveRelatedVar::NO);
          ASSERT_TRUE(dualProgram.has_value());
          for (const auto lexicographicReoptType : lexicographicReoptTypes) {
            IPOptStatistics<FloatingPointT> ipOptStatistics =
                runPrimalSimplexGomory<FloatingPointT, SimplexTraitsT>(
                    primalProgram, *dualProgram, simplexTableauType,
                    lexicographicReoptType, lpOptimizationType,
                    GomoryCutChoosingRule::FIRST);
            const auto gurobiLPOptStats =
                GurobiOptimizer("", modelFileMpsPath)
                    .optimize<FloatingPointT>(lpOptimizationType);
            this->compare(lpOptimizationType, lexicographicReoptType,
                          ipOptStatistics, gurobiLPOptStats,
                          IsDualProgramOptimized::YES);
            optStatsVec.push_back({ipOptStatistics, gurobiLPOptStats});
          }
        },
        AllBoundsMustBeSpecified::NO, DoVarsMustBeNonnegative::NO);
  }
};

TYPED_TEST_SUITE_P(PrimalSimplexGomoryTest);

TYPED_TEST_P(PrimalSimplexGomoryTest,
             runPrimalSimplexGomoryAndCompareWithGurobi) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_extended_statistics, true);
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/gomory_example_instances", 300,
      LPOptimizationType::LINEAR_RELAXATION, {LexicographicReoptType::MIN}));
}

TYPED_TEST_P(PrimalSimplexGomoryTest,
             runPrimalSimplexGomoryAndCompareWithGurobiSingleInstance) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_extended_statistics, true);
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/gomory_single_instance", 50,
      LPOptimizationType::LINEAR_RELAXATION, {LexicographicReoptType::MAX}));
}

TYPED_TEST_P(PrimalSimplexGomoryTest,
             runPrimalSimplexGomoryAndCompareWithGurobiWorkingInstancesINT) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_extended_statistics, true);
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/primal_cuts_working_instances", 500,
      LPOptimizationType::INTEGER_PROGRAM, {LexicographicReoptType::MAX}));
}

TYPED_TEST_P(PrimalSimplexGomoryTest,
             runPrimalSimplexGomoryAndCompareWithGurobiSingleInstanceINT) {
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_extended_statistics, true);
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/primal_cuts_single_instance", 1,
      LPOptimizationType::INTEGER_PROGRAM, {LexicographicReoptType::MAX}));
}

REGISTER_TYPED_TEST_SUITE_P(
    PrimalSimplexGomoryTest, runPrimalSimplexGomoryAndCompareWithGurobi,
    runPrimalSimplexGomoryAndCompareWithGurobiSingleInstance,
    runPrimalSimplexGomoryAndCompareWithGurobiWorkingInstancesINT,
    runPrimalSimplexGomoryAndCompareWithGurobiSingleInstanceINT);

using PrimalSimplexGomoryTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(PrimalSimplexGomoryTestSuite,
                               PrimalSimplexGomoryTest,
                               PrimalSimplexGomoryTypes);
