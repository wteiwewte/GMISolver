#include "src/Algorithms/DualSimplexGomory.h"
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
IPOptStatistics<T> runDualSimplexGomoryWithPrimalCuts(
    const LinearProgram<T> &linearProgram,
    const SimplexTableauType simplexTableauType,
    const LexicographicReoptType lexicographicReoptType,
    const LPOptimizationType lpOptimizationType,
    const GomoryCutChoosingRule gomoryCutChoosingRule) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  DualSimplexGomory<T, SimplexTraitsT> dualSimplexGomoryWithPrimalCuts(
      simplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option),
      absl::GetFlag(FLAGS_slack_cut_removal_condition), lexicographicReoptType);
  return dualSimplexGomoryWithPrimalCuts.run(lpOptimizationType,
                                             gomoryCutChoosingRule);
}

template <typename T>
class DualSimplexGomoryTest : public LPTestBase<T>, public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_simplex_tableau_types,
                  {SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE});
  }

  void
  testCase(const std::string &dualSimplexGomoryTestDirPath,
           const size_t basisSizeLimit,
           const LPOptimizationType lpOptimizationType,
           const std::vector<LexicographicReoptType> lexicographicReoptTypes) {
    using FloatingPointT = std::tuple_element_t<0, typename T::types>;
    using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
    using UnderlyingOptStatsT = IPOptStatistics<FloatingPointT>;
    this->template solveAndCompareInstancesFromSets<UnderlyingOptStatsT>(
        dualSimplexGomoryTestDirPath, basisSizeLimit, lpOptimizationType,
        [&](const auto &linearProgram,
            const SimplexTableauType simplexTableauType,
            const std::filesystem::path &modelFileMpsPath,
            std::vector<OptimizationStats<UnderlyingOptStatsT>> &optStatsVec) {
          for (const auto lexicographicReoptType : lexicographicReoptTypes) {
            IPOptStatistics<FloatingPointT> ipOptStatistics =
                runDualSimplexGomoryWithPrimalCuts<FloatingPointT,
                                                   SimplexTraitsT>(
                    linearProgram, simplexTableauType, lexicographicReoptType,
                    lpOptimizationType, GomoryCutChoosingRule::FIRST);
            const auto gurobiLPOptStats =
                GurobiOptimizer("", modelFileMpsPath)
                    .optimize<FloatingPointT>(lpOptimizationType);
            this->compare(lpOptimizationType, lexicographicReoptType,
                          ipOptStatistics, gurobiLPOptStats);
            optStatsVec.push_back({ipOptStatistics, gurobiLPOptStats});
          }
        });
  }
};

TYPED_TEST_SUITE_P(DualSimplexGomoryTest);

TYPED_TEST_P(DualSimplexGomoryTest,
             runDualSimplexWithLexReoptAndCompareWithGurobi) {
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/dual_simplex_working_instances", 200,
      LPOptimizationType::LINEAR_RELAXATION,
      {LexicographicReoptType::MIN, LexicographicReoptType::MAX}));
}
TYPED_TEST_P(DualSimplexGomoryTest, runDualSimplexGomoryAndCompareWithGurobi) {
  absl::SetFlag(&FLAGS_simplex_tableau_types, {SimplexTableauType::FULL});
  absl::SetFlag(&FLAGS_validate_simplex_option,
                ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
  absl::SetFlag(&FLAGS_extended_statistics, true);
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/gomory_example_instances", 50,
      LPOptimizationType::INTEGER_PROGRAM, {LexicographicReoptType::MAX}));
}

REGISTER_TYPED_TEST_SUITE_P(DualSimplexGomoryTest,
                            runDualSimplexGomoryAndCompareWithGurobi,
                            runDualSimplexWithLexReoptAndCompareWithGurobi);

using DualSimplexGomoryTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(DualSimplexGomoryTestSuite,
                               DualSimplexGomoryTest, DualSimplexGomoryTypes);
