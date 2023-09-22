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

template <typename T> struct PrimalSimplexOutput {
  LPOptStatistics<T> _phaseOneLpOptStats;
  std::optional<LPOptStatistics<T>> _phaseTwoLpOptStats;
};

template <typename T, typename SimplexTraitsT>
PrimalSimplexOutput<T>
runPrimalSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::PRIMAL,
      absl::GetFlag(FLAGS_use_product_form_of_inverse));
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

template <typename T> class PrimalSimplexTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TYPED_TEST_SUITE_P(PrimalSimplexTest);
TYPED_TEST_P(PrimalSimplexTest, runPrimalSimplexAndCompareWithGurobi) {
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

REGISTER_TYPED_TEST_SUITE_P(PrimalSimplexTest,
                            runPrimalSimplexAndCompareWithGurobi);

using PrimalSimplexTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
// using PrimalSimplexTypes = ::testing::Types<
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::SPARSE>>,
//     TypeTuple<double, SimplexTraits<double,
//     MatrixRepresentationType::NORMAL>>, TypeTuple<long double,
//               SimplexTraits<long double, MatrixRepresentationType::SPARSE>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(PrimalSimplexTestSuite, PrimalSimplexTest,
                               PrimalSimplexTypes);
