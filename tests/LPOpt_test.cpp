#include "Algorithms/SimplexTableau.h"
#include "Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"

#include <filesystem>

#include <absl/flags/flag.h>
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

ABSL_FLAG(
    int32_t, obj_value_logging_frequency, 0,
    "Current objective value should be logged every nth iteration of simplex");
ABSL_FLAG(int32_t, reinversion_frequency, 60,
          "Basis matrix should be reinverted every nth iteration of simplex");
ABSL_FLAG(bool, use_product_form_of_inverse, true,
          "Basis matrix inverse is represented via product form of inverse");

template <typename... Ts>
struct TypeTuple
{
  using types = std::tuple<Ts...>;
};

template <typename T, typename SimplexTraitsT>
LPOptStatistics<T> runDualSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, false, absl::GetFlag(FLAGS_use_product_form_of_inverse));
  return RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
                                   simplexTableau, DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
                                   absl::GetFlag(FLAGS_obj_value_logging_frequency),
                                   absl::GetFlag(FLAGS_reinversion_frequency))
                                   .run("");
}

template <typename T>
class LPOptTest : public ::testing::Test {
protected:
  void SetUp() override {

  }

};

TYPED_TEST_SUITE_P(LPOptTest);

TYPED_TEST_P(LPOptTest, runLPOptAndCompareWithGurobi) {
  using TypeTupleT = TypeParam;
  using FloatingPointT = std::tuple_element_t<0, typename TypeTupleT::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeTupleT::types>;

  for (const auto &lpModelFileEntry :
       std::filesystem::directory_iterator("../../tests/data")) //FIXME
  {
    auto linearProgram = MpsReader<FloatingPointT>::read(lpModelFileEntry.path());
    ASSERT_TRUE(linearProgram.has_value());

    const auto dualSimplexLpOptStats = runDualSimplexWithImplicitBounds<FloatingPointT, SimplexTraitsT>(*linearProgram);
    const auto gurobiLPOptStats = GurobiOptimizer("", lpModelFileEntry.path()).optimize(LPOptimizationType::LINEAR_RELAXATION);
    EXPECT_EQ(gurobiLPOptStats._optResult, dualSimplexLpOptStats._optResult);
    SPDLOG_INFO("MODEL {} GUROBI OPT {} SIMPLEX OPT {}", std::string{lpModelFileEntry.path().filename()}, gurobiLPOptStats._optimalValue,
                dualSimplexLpOptStats._optimalValue);
    EXPECT_NEAR(gurobiLPOptStats._optimalValue, dualSimplexLpOptStats._optimalValue, 0.00001);
  }
}

REGISTER_TYPED_TEST_SUITE_P(LPOptTest, runLPOptAndCompareWithGurobi);

typedef ::testing::Types<TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>,
                         TypeTuple<long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>,
                         TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>,
                         TypeTuple<long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>> FloatingTypes;
INSTANTIATE_TYPED_TEST_SUITE_P(LpOptSuite, LPOptTest, FloatingTypes);
