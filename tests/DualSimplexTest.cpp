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

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T>
runDualSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram,
                                 const SimplexTableauType simplexTableauType) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  return RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
             simplexTableau, reinversionManager,
             DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
             absl::GetFlag(FLAGS_obj_value_logging_frequency),
             absl::GetFlag(FLAGS_validate_simplex_option))
      .run("");
}

template <typename T>
class DualSimplexTest : public LPTestBase<T>, public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_simplex_tableau_types,
                  {SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE,
                   SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE});
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
          const SimplexTableauType simplexTableauType,
          const std::filesystem::path &modelFileMpsPath) {
        const auto dualSimplexLpOptStats =
            runDualSimplexWithImplicitBounds<FloatingPointT, SimplexTraitsT>(
                linearProgram, simplexTableauType);
        const auto gurobiLPOptStats =
            GurobiOptimizer("", modelFileMpsPath)
                .optimize<FloatingPointT>(
                    LPOptimizationType::LINEAR_RELAXATION);
        this->compareWithGurobi(dualSimplexLpOptStats, gurobiLPOptStats);
      });
}

REGISTER_TYPED_TEST_SUITE_P(DualSimplexTest,
                            runDualSimplexAndCompareWithGurobi);

using DualSimplexTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>,
    TypeTuple<long double,
              SimplexTraits<long double, MatrixRepresentationType::NORMAL>>>;
//    TypeTuple<double, SimplexTraits<double,
//    MatrixRepresentationType::SPARSE>>, TypeTuple<long double,
//              SimplexTraits<long double, MatrixRepresentationType::SPARSE>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(DualSimplexTestSuite, DualSimplexTest,
                               DualSimplexTypes);
