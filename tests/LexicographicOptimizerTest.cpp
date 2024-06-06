#include "Algorithms/LexicographicOptimizer.h"
#include "Algorithms/DualSimplex.h"
#include "Algorithms/DualSimplexPhaseOne.h"
#include "Algorithms/PrimalSimplex.h"
#include "Algorithms/PrimalSimplexPhaseOne.h"
#include "Algorithms/ReinversionManager.h"
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
    const SimplexTableauType simplexTableauType,
    const LexicographicReoptType lexicographicReoptType) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  DualSimplexPhaseOne<T, SimplexTraitsT> dualSimplexPhaseOne(
      simplexTableau, reinversionManager,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));

  auto phaseOneLpOptStats = dualSimplexPhaseOne.run();
  if (!phaseOneLpOptStats._phaseOneSucceeded) {
    SPDLOG_WARN("PHASE ONE OF {} ALGORITHM FAILED", dualSimplexPhaseOne.type());
    return {};
  }
  DualSimplex<T, SimplexTraitsT>(
      simplexTableau, reinversionManager,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option))
      .run("", DualPhase::TWO);

  return LexicographicOptimizer<T, SimplexTraitsT>(
             simplexTableau, reinversionManager,
             PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
             absl::GetFlag(FLAGS_obj_value_logging_frequency),
             absl::GetFlag(FLAGS_validate_simplex_option),
             lexicographicReoptType)
      .run("", SaveLexSolution::YES);
}

template <typename T, typename SimplexTraitsT>
LexReoptStatistics<T> runPrimalSimplexWithLexReopt(
    const LinearProgram<T> &linearProgram,
    const SimplexTableauType simplexTableauType,
    const LexicographicReoptType lexicographicReoptType) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::PRIMAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  PrimalSimplex<T, SimplexTraitsT> primalSimplex(
      simplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));

  PrimalSimplexPhaseOne<T, SimplexTraitsT> primalSimplexPhaseOne(
      simplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));
  auto phaseOneLpOptStats = primalSimplexPhaseOne.run();
  if (!phaseOneLpOptStats._phaseOneSucceeded) {
    SPDLOG_WARN("PHASE ONE OF {} ALGORITHM FAILED", primalSimplex.type());
    return {};
  }
  primalSimplex.runPhaseTwo();

  return LexicographicOptimizer<T, SimplexTraitsT>(
             simplexTableau, reinversionManager,
             PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
             absl::GetFlag(FLAGS_obj_value_logging_frequency),
             absl::GetFlag(FLAGS_validate_simplex_option),
             lexicographicReoptType)
      .run("", SaveLexSolution::YES);
}

template <typename T>
class LexicographicOptimizerTest : public LPTestBase<T>,
                                   public ::testing::Test {
protected:
  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex_option,
                  ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR);
    absl::SetFlag(&FLAGS_simplex_tableau_types,
                  {SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE,
                   SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE,
                   SimplexTableauType::FULL});
  }

  template <typename SimplexFunc>
  void testCase(const std::string &testDirPath, const size_t basisSizeLimit,
                SimplexFunc simplexFunc) {
    using FloatingPointT = std::tuple_element_t<0, typename T::types>;
    using UnderlyingOptStatsT = LexReoptStatistics<FloatingPointT>;
    const LPOptimizationType lpOptimizationType{
        LPOptimizationType::LINEAR_RELAXATION};
    this->template solveAndCompareInstancesFromSets<UnderlyingOptStatsT>(
        testDirPath, basisSizeLimit, lpOptimizationType,
        [&](const auto &linearProgram,
            const SimplexTableauType simplexTableauType,
            const std::filesystem::path &modelFileMpsPath,
            std::vector<OptimizationStats<UnderlyingOptStatsT>> &optStatsVec) {
          for (const auto lexicographicReoptType :
               {LexicographicReoptType::MIN, LexicographicReoptType::MAX}) {
            LexReoptStatistics<FloatingPointT> lexReoptStatistics = simplexFunc(
                linearProgram, simplexTableauType, lexicographicReoptType);

            GurobiOptimizer gurobiOptimizer("", modelFileMpsPath);
            gurobiOptimizer.prepareLexicographicObjectives(
                lexicographicReoptType);
            const auto gurobiLPOptStats =
                gurobiOptimizer.optimize<FloatingPointT>(lpOptimizationType);
            this->compareWithGurobi(
                lexicographicReoptType, lexReoptStatistics, gurobiLPOptStats,
                gurobiOptimizer.getSolutionVector<FloatingPointT>());
            optStatsVec.push_back({lexReoptStatistics, gurobiLPOptStats});
          }
        },
        false);
  }
};

TYPED_TEST_SUITE_P(LexicographicOptimizerTest);

TYPED_TEST_P(LexicographicOptimizerTest,
             runDualSimplexWithLexReoptAndCompareWithGurobi) {
  using FloatingPointT = std::tuple_element_t<0, typename TypeParam::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeParam::types>;
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/dual_simplex_working_instances", 1000,
      [](const auto &linearProgram, const SimplexTableauType simplexTableauType,
         const LexicographicReoptType lexicographicReoptType) {
        return runDualSimplexWithLexReopt<FloatingPointT, SimplexTraitsT>(
            linearProgram, simplexTableauType, lexicographicReoptType);
      }));
}
TYPED_TEST_P(LexicographicOptimizerTest,
             runPrimalSimplexWithLexReoptAndCompareWithGurobi) {
  using FloatingPointT = std::tuple_element_t<0, typename TypeParam::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename TypeParam::types>;
  absl::SetFlag(&FLAGS_extended_statistics, false);
  EXPECT_NO_FATAL_FAILURE(this->testCase(
      "../../tests/primal_simplex_working_instances", 500,
      [](const auto &linearProgram, const SimplexTableauType simplexTableauType,
         const LexicographicReoptType lexicographicReoptType) {
        return runPrimalSimplexWithLexReopt<FloatingPointT, SimplexTraitsT>(
            linearProgram, simplexTableauType, lexicographicReoptType);
      }));
}

REGISTER_TYPED_TEST_SUITE_P(LexicographicOptimizerTest,
                            runDualSimplexWithLexReoptAndCompareWithGurobi,
                            runPrimalSimplexWithLexReoptAndCompareWithGurobi);

using LexicographicOptimizerTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(LexicographicOptimizerTestSuite,
                               LexicographicOptimizerTest,
                               LexicographicOptimizerTypes);
