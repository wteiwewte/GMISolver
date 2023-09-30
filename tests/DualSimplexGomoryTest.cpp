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
    const LexicographicReoptType lexicographicReoptType,
    const LPOptimizationType lpOptimizationType,
    const GomoryCutChoosingRule gomoryCutChoosingRule) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL,
      absl::GetFlag(FLAGS_use_product_form_of_inverse));
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  DualSimplexGomory<T, SimplexTraitsT> dualSimplexGomoryWithPrimalCuts(
      simplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));
  return dualSimplexGomoryWithPrimalCuts.run(
      lexicographicReoptType, lpOptimizationType, gomoryCutChoosingRule);
}

template <typename T>
class DualSimplexGomoryTest : public LPTestBase<T>, public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_use_product_form_of_inverse, true);
  }
};

TYPED_TEST_SUITE_P(DualSimplexGomoryTest);

TYPED_TEST_P(DualSimplexGomoryTest,
             runDualSimplexWithLexReoptAndCompareWithGurobi) {
  using FloatingPointT = std::tuple_element_t<0, typename TypeParam::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeParam::types>;
  constexpr auto DUAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/dual_simplex_working_instances";
  constexpr size_t DUAL_SIMPLEX_BASIS_SIZE_LIMIT = 1000;
  const LPOptimizationType lpOptimizationType{
      LPOptimizationType::LINEAR_RELAXATION};
  this->solveAndCompareInstancesFromSets(
      DUAL_SIMPLEX_TEST_DIR_PATH, DUAL_SIMPLEX_BASIS_SIZE_LIMIT,
      lpOptimizationType,
      [&](const auto &linearProgram,
          const std::filesystem::path &modelFileMpsPath) {
        for (const auto lexicographicReoptType :
             {LexicographicReoptType::MIN, LexicographicReoptType::MAX}) {
          IPOptStatistics<FloatingPointT> ipOptStatistics =
              runDualSimplexGomoryWithPrimalCuts<FloatingPointT,
                                                 SimplexTraitsT>(
                  linearProgram, lexicographicReoptType, lpOptimizationType,
                  GomoryCutChoosingRule::FIRST);
          const auto gurobiLPOptStats =
              GurobiOptimizer("", modelFileMpsPath)
                  .optimize<FloatingPointT>(lpOptimizationType);
          this->compare(lpOptimizationType, lexicographicReoptType,
                        ipOptStatistics, gurobiLPOptStats);
        }
      });
}
TYPED_TEST_P(DualSimplexGomoryTest, runDualSimplexGomoryAndCompareWithGurobi) {
  using FloatingPointT = std::tuple_element_t<0, typename TypeParam::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeParam::types>;
  constexpr auto DUAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/gomory_example_instances";
  //      "../../tests/dual_simplex_working_instances";

  constexpr size_t DUAL_SIMPLEX_BASIS_SIZE_LIMIT = 25;
  const LPOptimizationType lpOptimizationType{
      LPOptimizationType::INTEGER_PROGRAM};
  this->solveAndCompareInstancesFromSets(
      DUAL_SIMPLEX_TEST_DIR_PATH, DUAL_SIMPLEX_BASIS_SIZE_LIMIT,
      lpOptimizationType,
      [&](const auto &linearProgram,
          const std::filesystem::path &modelFileMpsPath) {
        for (const auto lexicographicReoptType :
             {LexicographicReoptType::MAX}) {
          IPOptStatistics<FloatingPointT> ipOptStatistics =
              runDualSimplexGomoryWithPrimalCuts<FloatingPointT,
                                                 SimplexTraitsT>(
                  linearProgram, lexicographicReoptType, lpOptimizationType,
                  GomoryCutChoosingRule::FIRST);
          const auto gurobiLPOptStats =
              GurobiOptimizer("", modelFileMpsPath)
                  .optimize<FloatingPointT>(lpOptimizationType);
          this->compare(lpOptimizationType, lexicographicReoptType,
                        ipOptStatistics, gurobiLPOptStats);
        }
      });
}

REGISTER_TYPED_TEST_SUITE_P(DualSimplexGomoryTest,
                            runDualSimplexGomoryAndCompareWithGurobi,
                            runDualSimplexWithLexReoptAndCompareWithGurobi);

using DualSimplexGomoryTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(DualSimplexTestSuite, DualSimplexGomoryTest,
                               DualSimplexGomoryTypes);
