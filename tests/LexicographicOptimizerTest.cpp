#include "Algorithms/LexicographicOptimizer.h"
#include "Algorithms/ReinversionManager.h"
#include "Algorithms/RevisedDualSimplexPFIBounds.h"
#include "Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"
#include "tests/CommonDefs.h"
#include "tests/LPTestBase.h"

#include <filesystem>

#include <gtest/gtest.h>

template <typename T, typename SimplexTraitsT>
LexReoptStatistics<T> runDualSimplexWithLexReopt(
    const LinearProgram<T> &linearProgram,
    const LexicographicReoptType lexicographicReoptType) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL,
      absl::GetFlag(FLAGS_use_product_form_of_inverse));
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
      simplexTableau, reinversionManager,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option))
      .run("");
  return LexicographicOptimizer<T, SimplexTraitsT>(
             simplexTableau, reinversionManager,
             PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
             absl::GetFlag(FLAGS_obj_value_logging_frequency),
             absl::GetFlag(FLAGS_validate_simplex_option))
      .run(lexicographicReoptType, "", true);
}

template <typename T>
class LexicographicOptimizerTest : public LPTestBase<T>,
                                   public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_use_product_form_of_inverse, true);
  }
};

TYPED_TEST_SUITE_P(LexicographicOptimizerTest);

TYPED_TEST_P(LexicographicOptimizerTest,
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
          LexReoptStatistics<FloatingPointT> lexReoptStatistics =
              runDualSimplexWithLexReopt<FloatingPointT, SimplexTraitsT>(
                  linearProgram, lexicographicReoptType);

          GurobiOptimizer gurobiOptimizer("", modelFileMpsPath);
          gurobiOptimizer.prepareLexicographicObjectives(
              lexicographicReoptType);
          const auto gurobiLPOptStats =
              gurobiOptimizer.optimize<FloatingPointT>(lpOptimizationType);

          const auto gurobiSolution =
              gurobiOptimizer.getSolutionVector<FloatingPointT>();
          const auto &lexOptimizerSolution = lexReoptStatistics._solution;
          ASSERT_EQ(gurobiSolution.size(), lexOptimizerSolution.size());
          for (int varIdx = 0; varIdx < gurobiSolution.size(); ++varIdx) {
            EXPECT_NEAR(gurobiSolution[varIdx], lexOptimizerSolution[varIdx],
                        0.00001);
          }
          this->compareWithGurobi(lexicographicReoptType, lexReoptStatistics,
                                  gurobiLPOptStats);
        }
      });
}

REGISTER_TYPED_TEST_SUITE_P(LexicographicOptimizerTest,
                            runDualSimplexWithLexReoptAndCompareWithGurobi);

// using LexicographicOptimizerTypes = ::testing::Types<
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>>;
using LexicographicOptimizerTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>,
    TypeTuple<long double,
              SimplexTraits<long double, MatrixRepresentationType::NORMAL>>>;
//    TypeTuple<double, SimplexTraits<double,
//    MatrixRepresentationType::SPARSE>>, TypeTuple<long double,
//              SimplexTraits<long double, MatrixRepresentationType::SPARSE>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(LexicographicOptimizerTestSuite,
                               LexicographicOptimizerTest,
                               LexicographicOptimizerTypes);
