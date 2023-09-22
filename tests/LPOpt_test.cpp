#include "Algorithms/RevisedDualSimplexPFIBounds.h"
#include "Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "Algorithms/SimplexTableau.h"
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
ABSL_FLAG(bool, validate_simplex, false, "Validate simplex implementations");

template <typename... Ts> struct TypeTuple {
  using types = std::tuple<Ts...>;
};

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

template <typename T> struct PrimalSimplexOutput {
  LPOptStatistics<T> _phaseOneLpOptStats;
  std::optional<LPOptStatistics<T>> _phaseTwoLpOptStats;
};

template <typename T, typename SimplexTraitsT>
PrimalSimplexOutput<T>
runPrimalSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, true, absl::GetFlag(FLAGS_use_product_form_of_inverse));
  RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>
      revisedPrimalSimplexPfiBounds(
          simplexTableau,
          PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
          absl::GetFlag(FLAGS_obj_value_logging_frequency),
          absl::GetFlag(FLAGS_reinversion_frequency),
          absl::GetFlag(FLAGS_validate_simplex));
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

template <typename T> class LPOptTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TYPED_TEST_SUITE_P(LPOptTest);

TYPED_TEST_P(LPOptTest, runDualSimplexAndCompareWithGurobi) {
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
TYPED_TEST_P(LPOptTest, runPrimalSimplexAndCompareWithGurobi) {
  constexpr auto PRIMAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/primal_simplex_working_instances";
  constexpr size_t PRIMAL_SIMPLEX_BASIS_SIZE_LIMIT = 200;
  using TypeTupleT = TypeParam;
  using FloatingPointT = std::tuple_element_t<0, typename TypeTupleT::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeTupleT::types>;

  for (const auto &lpModelSetDirectory :
       std::filesystem::directory_iterator(PRIMAL_SIMPLEX_TEST_DIR_PATH)) {
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

      if (linearProgram->getRowInfos().size() > PRIMAL_SIMPLEX_BASIS_SIZE_LIMIT)
        continue;

      const auto primalSimplexOutput =
          runPrimalSimplexWithImplicitBounds<FloatingPointT, SimplexTraitsT>(
              *linearProgram);
      const auto gurobiLPOptStats =
          GurobiOptimizer("", lpModelFileEntry.path())
              .optimize(LPOptimizationType::LINEAR_RELAXATION);
      if (primalSimplexOutput._phaseOneLpOptStats._optResult ==
          LPOptimizationResult::INFEASIBLE) {
        const std::set<LPOptimizationResult> infeasibleResults{
            LPOptimizationResult::INFEASIBLE,
            LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED};
        ASSERT_TRUE(infeasibleResults.contains(gurobiLPOptStats._optResult));
      } else {
        ASSERT_TRUE(primalSimplexOutput._phaseTwoLpOptStats.has_value());
        const auto &phaseTwoLpOptStats =
            *primalSimplexOutput._phaseTwoLpOptStats;
        ASSERT_EQ(gurobiLPOptStats._optResult, phaseTwoLpOptStats._optResult);
        SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
        SPDLOG_INFO("SIMPLEX OPT {}", phaseTwoLpOptStats._optimalValue);
        EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                    phaseTwoLpOptStats._optimalValue, 0.00001);
      }
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(LPOptTest, runDualSimplexAndCompareWithGurobi,
                            runPrimalSimplexAndCompareWithGurobi);

using SimplexTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
// using SimplexTypes = ::testing::Types<
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::SPARSE>>,
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::SPARSE>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(LpOptSuite, LPOptTest, SimplexTypes);
