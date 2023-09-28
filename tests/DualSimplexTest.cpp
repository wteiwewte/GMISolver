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

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T>
runDualSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL,
      absl::GetFlag(FLAGS_use_product_form_of_inverse));
  return RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
             simplexTableau, DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
             absl::GetFlag(FLAGS_obj_value_logging_frequency),
             absl::GetFlag(FLAGS_reinversion_frequency),
             absl::GetFlag(FLAGS_validate_simplex_option))
      .run("");
}

template <typename T>
class DualSimplexTest : public LPTestBase<T>, public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_use_product_form_of_inverse, true);
  }
};

TYPED_TEST_SUITE_P(DualSimplexTest);

TYPED_TEST_P(DualSimplexTest, runDualSimplexAndCompareWithGurobi) {
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
        const auto dualSimplexLpOptStats =
            runDualSimplexWithImplicitBounds<FloatingPointT, SimplexTraitsT>(
                linearProgram);
        const auto gurobiLPOptStats =
            GurobiOptimizer("", modelFileMpsPath)
                .optimize(LPOptimizationType::LINEAR_RELAXATION);
        this->compareWithGurobi(dualSimplexLpOptStats, gurobiLPOptStats);
      });
}

REGISTER_TYPED_TEST_SUITE_P(DualSimplexTest,
                            runDualSimplexAndCompareWithGurobi);

using DualSimplexTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
// using DualSimplexTypes = ::testing::Types<
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::SPARSE>>,
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::SPARSE>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(DualSimplexTestSuite, DualSimplexTest,
                               DualSimplexTypes);
