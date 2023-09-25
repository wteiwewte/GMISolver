#include "Algorithms/DualSimplexGomory.h"
#include "Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"
#include "tests/CommonDefs.h"

#include <filesystem>

#include <gtest/gtest.h>

template <typename T, typename SimplexTraitsT>
IPOptStatistics<T>
runDualSimplexGomoryWithPrimalCuts(const LinearProgram<T> &linearProgram) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL,
      absl::GetFlag(FLAGS_use_product_form_of_inverse));
  DualSimplexGomory<T, SimplexTraitsT> dualSimplexGomoryWithPrimalCuts(
      simplexTableau,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_reinversion_frequency),
      absl::GetFlag(FLAGS_validate_simplex));
  return dualSimplexGomoryWithPrimalCuts.run();
}

template <typename T> class DualSimplexGomoryTest : public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex, ValidateSimplex::YES);
    absl::SetFlag(&FLAGS_use_product_form_of_inverse, true);
  }
};

TYPED_TEST_SUITE_P(DualSimplexGomoryTest);

TYPED_TEST_P(DualSimplexGomoryTest, runDualSimplexGomoryAndCompareWithGurobi) {
  constexpr auto DUAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/dual_simplex_working_instances";
  //  constexpr auto DUAL_SIMPLEX_NOT_WORKING_TEST_DIR_PATH =
  //      "../../tests/dual_simplex_not_working_instances/";
  constexpr size_t DUAL_SIMPLEX_BASIS_SIZE_LIMIT = 1000;

  using TypeTupleT = TypeParam;
  using FloatingPointT = std::tuple_element_t<0, typename TypeTupleT::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeTupleT::types>;

  int processedModelsCount = 0;
  int tooBigModelsCount = 0;
  int modelsWithNoBoundsCount = 0;
  int modelsWithNoAllIntegerVariables = 0;
  int totalModelsCount = 0;
  for (const auto &lpModelSetDirectory :
       std::filesystem::directory_iterator(DUAL_SIMPLEX_TEST_DIR_PATH)) {
    if (!lpModelSetDirectory.is_directory())
      continue;

    SPDLOG_INFO("MODEL SET DIRECTORY {}",
                std::string{lpModelSetDirectory.path().filename()});
    for (const auto &lpModelFileEntry :
         std::filesystem::directory_iterator(lpModelSetDirectory)) {
      ++totalModelsCount;
      SPDLOG_INFO("MODEL {}", std::string{lpModelFileEntry.path().filename()});
      auto linearProgram =
          MpsReader<FloatingPointT>::read(lpModelFileEntry.path());
      ASSERT_TRUE(linearProgram.has_value());

      if (linearProgram->getRowInfos().size() > DUAL_SIMPLEX_BASIS_SIZE_LIMIT ||
          linearProgram->getVariableInfos().size() >
              2 * DUAL_SIMPLEX_BASIS_SIZE_LIMIT) {
        ++tooBigModelsCount;
        continue;
      }

      if (!linearProgram->checkIfAllBoundsAreSpeficied()) {
        ++modelsWithNoBoundsCount;
        SPDLOG_INFO(
            "SKIPPING MODEL {} BECAUSE NOT ALL VARIABLES HAVE BOUNDS SPECIFIED",
            std::string{lpModelFileEntry.path().filename()});
        continue;
      }

      if (!linearProgram->isPureIP()) {
        ++modelsWithNoAllIntegerVariables;
        SPDLOG_INFO("SKIPPING MODEL {} BECAUSE NOT ALL VARIABLES ARE INTEGER",
                    std::string{lpModelFileEntry.path().filename()});
        continue;
      }

      ++processedModelsCount;

      IPOptStatistics<FloatingPointT> ipOptStatistics =
          runDualSimplexGomoryWithPrimalCuts<FloatingPointT, SimplexTraitsT>(
              *linearProgram);
      const auto gurobiLPOptStats =
          GurobiOptimizer("", lpModelFileEntry.path())
              .optimize(LPOptimizationType::LINEAR_RELAXATION);
      const auto &firstRelaxationOptStats =
          ipOptStatistics._lpRelaxationStats.front();
      const auto &firstRelaxationLPOptStats =
          firstRelaxationOptStats._relaxationOptStats;
      ASSERT_EQ(gurobiLPOptStats._optResult,
                firstRelaxationLPOptStats._optResult);
      if (firstRelaxationLPOptStats._optResult ==
          LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
        SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
        SPDLOG_INFO("SIMPLEX OPT AFTER INITIAL OPT {}",
                    firstRelaxationLPOptStats._optimalValue);
        EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                    firstRelaxationLPOptStats._optimalValue, 0.00001);

        const auto &optimalValueAfterFirstLexReopt =
            firstRelaxationOptStats._lexicographicReoptStats
                ._objectiveValueAfterLexReopt;
        SPDLOG_INFO("SIMPLEX OPT AFTER LEXICOGRAPHIC REOPTS {}",
                    optimalValueAfterFirstLexReopt);
        EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                    optimalValueAfterFirstLexReopt, 0.00001);
      }
    }
  }
  SPDLOG_INFO("PROCESSED {} MODELS FROM {} TOTAL ({} TOO BIG, {} WITHOUT "
              "BOUNDS, {} NOT PURE IP)",
              processedModelsCount, totalModelsCount, tooBigModelsCount,
              modelsWithNoBoundsCount, modelsWithNoAllIntegerVariables);
}

REGISTER_TYPED_TEST_SUITE_P(DualSimplexGomoryTest,
                            runDualSimplexGomoryAndCompareWithGurobi);

using DualSimplexGomoryTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(DualSimplexTestSuite, DualSimplexGomoryTest,
                               DualSimplexGomoryTypes);
