#include "Algorithms/RevisedDualSimplexPFIBounds.h"
#include "Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"
#include "tests/CommonDefs.h"

#include <filesystem>

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T>
runDualSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, false, absl::GetFlag(FLAGS_use_product_form_of_inverse));
  return RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
             simplexTableau, DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
             absl::GetFlag(FLAGS_obj_value_logging_frequency),
             absl::GetFlag(FLAGS_reinversion_frequency),
             absl::GetFlag(FLAGS_validate_simplex))
      .run("");
}

template <typename T> class DualSimplexTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TYPED_TEST_SUITE_P(DualSimplexTest);

TYPED_TEST_P(DualSimplexTest, runDualSimplexAndCompareWithGurobi) {
  constexpr auto DUAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/dual_simplex_working_instances";
  //  constexpr auto DUAL_SIMPLEX_NOT_WORKING_TEST_DIR_PATH =
  //      "../../tests/dual_simplex_not_working_instances/";
  constexpr size_t DUAL_SIMPLEX_BASIS_SIZE_LIMIT = 500;

  using TypeTupleT = TypeParam;
  using FloatingPointT = std::tuple_element_t<0, typename TypeTupleT::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeTupleT::types>;

  for (const auto &lpModelSetDirectory :
       std::filesystem::directory_iterator(DUAL_SIMPLEX_TEST_DIR_PATH)) {
    if (!lpModelSetDirectory.is_directory())
      continue;

    SPDLOG_INFO("MODEL SET DIRECTORY {}",
                std::string{lpModelSetDirectory.path().filename()});
    for (const auto &lpModelFileEntry :
         std::filesystem::directory_iterator(lpModelSetDirectory)) {
      SPDLOG_INFO("MODEL {}", std::string{lpModelFileEntry.path().filename()});
      auto linearProgram =
          MpsReader<FloatingPointT>::read(lpModelFileEntry.path());
      ASSERT_TRUE(linearProgram.has_value());

      if (linearProgram->getRowInfos().size() > DUAL_SIMPLEX_BASIS_SIZE_LIMIT)
        continue;

      const auto dualSimplexLpOptStats =
          runDualSimplexWithImplicitBounds<FloatingPointT, SimplexTraitsT>(
              *linearProgram);
      const auto gurobiLPOptStats =
          GurobiOptimizer("", lpModelFileEntry.path())
              .optimize(LPOptimizationType::LINEAR_RELAXATION);
      ASSERT_EQ(gurobiLPOptStats._optResult, dualSimplexLpOptStats._optResult);
      if (dualSimplexLpOptStats._optResult ==
          LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
        SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
        SPDLOG_INFO("SIMPLEX OPT {}", dualSimplexLpOptStats._optimalValue);
        EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                    dualSimplexLpOptStats._optimalValue, 0.00001);
      }
    }
  }
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
