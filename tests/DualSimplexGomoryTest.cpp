#include "Algorithms/DualSimplexGomory.h"
#include "Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/MpsReader.h"
#include "tests/CommonDefs.h"

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
  DualSimplexGomory<T, SimplexTraitsT> dualSimplexGomoryWithPrimalCuts(
      simplexTableau,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_reinversion_frequency),
      absl::GetFlag(FLAGS_validate_simplex));
  return dualSimplexGomoryWithPrimalCuts.run(
      lexicographicReoptType, lpOptimizationType, gomoryCutChoosingRule);
}

struct InstanceSetStats {
  int _processedModelsCount = 0;
  int _tooBigModelsCount = 0;
  int _modelsWithNoBoundsCount = 0;
  int _modelsWithNoAllIntegerVariables = 0;
  int _modelsWithNoAllIntegerCoeffs = 0;
  int _modelsWithNoAllNonnegativeVariables = 0;
  int _totalModelsCount = 0;
};

template <typename T> class DualSimplexGomoryTest : public ::testing::Test {
protected:
  using FloatingPointT = std::tuple_element_t<0, typename T::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  void SetUp() override {
    absl::SetFlag(&FLAGS_validate_simplex, ValidateSimplex::YES);
    absl::SetFlag(&FLAGS_use_product_form_of_inverse, true);
  }

  bool isInstanceSuitable(const size_t basisSizeLimit,
                          const LinearProgram<FloatingPointT> &linearProgram,
                          const std::string &modelName,
                          const bool isRelaxationOptType,
                          InstanceSetStats &instanceSetStats) const {
    if (linearProgram.getRowInfos().size() > basisSizeLimit ||
        linearProgram.getVariableInfos().size() > 2 * basisSizeLimit) {
      ++instanceSetStats._tooBigModelsCount;
      return false;
    }

    if (!linearProgram.checkIfAllBoundsAreSpeficied()) {
      ++instanceSetStats._modelsWithNoBoundsCount;
      SPDLOG_INFO("SKIPPING MODEL {} BECAUSE NOT ALL VARIABLES HAVE BOUNDS "
                  "SPECIFIED",
                  modelName);
      return false;
    }

    if (!isRelaxationOptType) {
      if (!linearProgram.isPureIP()) {
        ++instanceSetStats._modelsWithNoAllIntegerVariables;
        SPDLOG_INFO("SKIPPING MODEL {} BECAUSE NOT ALL VARIABLES ARE INTEGER",
                    modelName);
        return false;
      }
      if (!linearProgram.allCoefficientsAreIntegers()) {
        ++instanceSetStats._modelsWithNoAllIntegerCoeffs;
        SPDLOG_INFO(
            "SKIPPING MODEL {} BECAUSE NOT ALL COEFFICIENTS ARE INTEGER",
            modelName);
        return false;
      }

      if (!linearProgram.allVariablesAreNonnegative()) {
        ++instanceSetStats._modelsWithNoAllNonnegativeVariables;
        SPDLOG_INFO(
            "SKIPPING MODEL {} BECAUSE NOT ALL VARIABLES ARE NONNEGATIVE",
            modelName);
        return false;
      }
    }

    ++instanceSetStats._processedModelsCount;
    return true;
  }

  void solveAndCompareInstancesFromSets(
      const std::string &dirPath, const size_t basisSizeLimit,
      const LPOptimizationType lpOptimizationType,
      const std::vector<LexicographicReoptType> lexReoptTypes) {
    InstanceSetStats instanceSetStats;
    const bool isRelaxationOptType =
        lpOptimizationType == LPOptimizationType::LINEAR_RELAXATION;
    for (const auto &lpModelSetDirectory :
         std::filesystem::directory_iterator(dirPath)) {
      if (!lpModelSetDirectory.is_directory())
        continue;

      SPDLOG_INFO("MODEL SET DIRECTORY {}",
                  std::string{lpModelSetDirectory.path().filename()});
      for (const auto &lpModelFileEntry :
           std::filesystem::directory_iterator(lpModelSetDirectory)) {
        ++instanceSetStats._totalModelsCount;
        const std::string modelName = lpModelFileEntry.path().filename();
        SPDLOG_INFO("MODEL {}", modelName);
        auto linearProgram =
            MpsReader<FloatingPointT>::read(lpModelFileEntry.path());
        ASSERT_TRUE(linearProgram.has_value());

        if (isInstanceSuitable(basisSizeLimit, *linearProgram, modelName,
                               isRelaxationOptType, instanceSetStats)) {
          for (const auto lexicographicReoptType : lexReoptTypes) {
            IPOptStatistics<FloatingPointT> ipOptStatistics =
                runDualSimplexGomoryWithPrimalCuts<FloatingPointT,
                                                   SimplexTraitsT>(
                    *linearProgram, lexicographicReoptType, lpOptimizationType,
                    GomoryCutChoosingRule::FIRST);
            const auto gurobiLPOptStats =
                GurobiOptimizer("", lpModelFileEntry.path())
                    .optimize(lpOptimizationType);
            compare(lpOptimizationType, lexicographicReoptType, ipOptStatistics,
                    gurobiLPOptStats);
          }
        }
      }
    }
    SPDLOG_INFO("PROCESSED {} MODELS FROM {} TOTAL",
                instanceSetStats._processedModelsCount,
                instanceSetStats._totalModelsCount);
    SPDLOG_INFO("{} TOO BIG", instanceSetStats._tooBigModelsCount);
    SPDLOG_INFO("{} WITHOUT BOUNDS", instanceSetStats._modelsWithNoBoundsCount);
    if (isRelaxationOptType) {
      ASSERT_EQ(instanceSetStats._totalModelsCount,
                instanceSetStats._processedModelsCount +
                    instanceSetStats._tooBigModelsCount +
                    instanceSetStats._modelsWithNoBoundsCount);
    } else {
      SPDLOG_INFO("{} NOT PURE IP",
                  instanceSetStats._modelsWithNoAllIntegerVariables);
      SPDLOG_INFO("{} NOT ALL COEFFS INTEGER",
                  instanceSetStats._modelsWithNoAllIntegerCoeffs);
      SPDLOG_INFO("{} NOT ALL VARIABLES ARE NONNEGATIVE",
                  instanceSetStats._modelsWithNoAllNonnegativeVariables);
      ASSERT_EQ(instanceSetStats._totalModelsCount,
                instanceSetStats._processedModelsCount +
                    instanceSetStats._tooBigModelsCount +
                    instanceSetStats._modelsWithNoBoundsCount +
                    instanceSetStats._modelsWithNoAllIntegerVariables +
                    instanceSetStats._modelsWithNoAllIntegerCoeffs +
                    instanceSetStats._modelsWithNoAllNonnegativeVariables);
    }
  }

  void compare(const LPOptimizationType lpOptimizationType,
               const LexicographicReoptType lexicographicReoptType,
               const IPOptStatistics<FloatingPointT> &ipOptStatistics,
               const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    if (lpOptimizationType == LPOptimizationType::LINEAR_RELAXATION) {
      compareLPRelaxation(lexicographicReoptType, ipOptStatistics,
                          gurobiLPOptStats);
    } else {
      compareIP(lexicographicReoptType, ipOptStatistics, gurobiLPOptStats);
    }
  }

  void
  compareLPRelaxation(const LexicographicReoptType lexicographicReoptType,
                      const IPOptStatistics<FloatingPointT> &ipOptStatistics,
                      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    const auto &firstRelaxationOptStats =
        ipOptStatistics._lpRelaxationStats.front();
    compareWithGurobi(lexicographicReoptType, firstRelaxationOptStats,
                      gurobiLPOptStats);
  }
  void compareIP(const LexicographicReoptType lexicographicReoptType,
                 const IPOptStatistics<FloatingPointT> &ipOptStatistics,
                 const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    const auto &lastRelaxationOptStats =
        ipOptStatistics._lpRelaxationStats.back();
    compareWithGurobi(lexicographicReoptType, lastRelaxationOptStats,
                      gurobiLPOptStats);
    checkIfSolutionIsInteger(ipOptStatistics);
  }

  void checkIfSolutionIsInteger(
      const IPOptStatistics<FloatingPointT> &ipOptStatistics) {
    const auto &optSolution = ipOptStatistics._optimalSolution;
    EXPECT_TRUE(
        std::all_of(optSolution.begin(), optSolution.end(), [](const auto &x) {
          return NumericalTraitsT ::isInteger(x);
        }));
  }

  void compareWithGurobi(
      const LexicographicReoptType lexicographicReoptType,
      const LPRelaxationStatistics<FloatingPointT> &lpRelaxationStatistics,
      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    const auto &relaxationOptStats = lpRelaxationStatistics._relaxationOptStats;
    ASSERT_EQ(gurobiLPOptStats._optResult, relaxationOptStats._optResult);
    if (relaxationOptStats._optResult ==
        LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
      SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
      SPDLOG_INFO("SIMPLEX OPT AFTER INITIAL OPT {}",
                  relaxationOptStats._optimalValue);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  relaxationOptStats._optimalValue,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);

      const auto &optimalValueAfterFirstLexReopt =
          lpRelaxationStatistics._lexicographicReoptStats
              ._objectiveValueAfterLexReopt;
      SPDLOG_INFO("SIMPLEX OPT AFTER LEXICOGRAPHIC {} REOPTS {}",
                  lexicographicReoptTypeToStr(lexicographicReoptType),
                  optimalValueAfterFirstLexReopt);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  optimalValueAfterFirstLexReopt,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }
};

TYPED_TEST_SUITE_P(DualSimplexGomoryTest);

TYPED_TEST_P(DualSimplexGomoryTest,
             runDualSimplexWithLexReoptAndCompareWithGurobi) {
  constexpr auto DUAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/dual_simplex_working_instances";
  constexpr size_t DUAL_SIMPLEX_BASIS_SIZE_LIMIT = 1000;
  this->solveAndCompareInstancesFromSets(
      DUAL_SIMPLEX_TEST_DIR_PATH, DUAL_SIMPLEX_BASIS_SIZE_LIMIT,
      LPOptimizationType::LINEAR_RELAXATION,
      {LexicographicReoptType::MIN, LexicographicReoptType::MAX});
}
TYPED_TEST_P(DualSimplexGomoryTest, runDualSimplexGomoryAndCompareWithGurobi) {
  constexpr auto DUAL_SIMPLEX_TEST_DIR_PATH =
      "../../tests/gomory_example_instances";
  //      "../../tests/dual_simplex_working_instances";

  constexpr size_t DUAL_SIMPLEX_BASIS_SIZE_LIMIT = 25;
  this->solveAndCompareInstancesFromSets(
      DUAL_SIMPLEX_TEST_DIR_PATH, DUAL_SIMPLEX_BASIS_SIZE_LIMIT,
      LPOptimizationType::INTEGER_PROGRAM, {LexicographicReoptType::MIN});
}

REGISTER_TYPED_TEST_SUITE_P(DualSimplexGomoryTest,
                            runDualSimplexGomoryAndCompareWithGurobi,
                            runDualSimplexWithLexReoptAndCompareWithGurobi);

using DualSimplexGomoryTypes = ::testing::Types<
    TypeTuple<double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(DualSimplexTestSuite, DualSimplexGomoryTest,
                               DualSimplexGomoryTypes);
