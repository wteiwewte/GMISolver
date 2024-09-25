#ifndef GMISOLVER_LPTESTBASE_H
#define GMISOLVER_LPTESTBASE_H

#include "src/Algorithms/SimplexTableau.h"
#include "src/DataModel/EnumTypes.h"
#include "src/Util/IPOptStatistics.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/OptStatisticsPrinter.h"
#include "src/Util/SpdlogHeader.h"
#include "tests/CommonDefs.h"

#include <tuple>

#include <gtest/gtest.h>

struct InstanceSetStats {
  int _processedModelsCount = 0;
  int _tooBigModelsCount = 0;
  int _modelsWithNoBoundsCount = 0;
  int _modelsWithNotAllIntegerVariables = 0;
  int _modelsWithNoAllIntegerCoeffs = 0;
  int _modelsWithNoAllNonnegativeVariables = 0;
  int _totalModelsCount = 0;
};

template <typename UnderlyingStatsT> struct OptimizationStats {
  using T = UnderlyingStatsT::Type;
  UnderlyingStatsT _optimizationStats;
  LPOptStatistics<T> _gurobiStats;
};

template <typename T> struct LPTestBase {
  using FloatingPointT = std::tuple_element_t<0, typename T::types>;
  using SimplexTraitsT = std::tuple_element_t<1, typename T::types>;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  bool isInstanceSuitable(const size_t basisSizeLimit,
                          const LinearProgram<FloatingPointT> &linearProgram,
                          const std::string &modelName,
                          const bool isRelaxationOptType,
                          const bool allBoundsMustBeSpecified,
                          InstanceSetStats &instanceSetStats) const {
    //        if (modelName != "air04.mps") {
    //          return false;
    //        }
    ++instanceSetStats._totalModelsCount;
    if (linearProgram.getRowInfos().size() > basisSizeLimit ||
        linearProgram.getVariableInfos().size() > 2 * basisSizeLimit) {
      ++instanceSetStats._tooBigModelsCount;
      SPDLOG_INFO("SKIPPING MODEL {} BECAUSE OF TOO BIG SIZE", modelName);
      return false;
    }

    if (allBoundsMustBeSpecified &&
        !linearProgram.checkIfAllBoundsAreSpeficied()) {
      ++instanceSetStats._modelsWithNoBoundsCount;
      SPDLOG_INFO("SKIPPING MODEL {} BECAUSE NOT ALL VARIABLES HAVE BOUNDS "
                  "SPECIFIED",
                  modelName);
      return false;
    }

    if (!isRelaxationOptType) {
      if (!linearProgram.isPureIP()) {
        ++instanceSetStats._modelsWithNotAllIntegerVariables;
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

  template <typename UnderlyingStatsT, typename LPOptFunc>
  void solveAndCompareInstancesFromSets(
      const std::string &dirPath, const size_t basisSizeLimit,
      const LPOptimizationType lpOptimizationType, LPOptFunc lpOptFunc,
      const bool allBoundsMustBeSpecified = true) {
    using OptimizationStatsT = OptimizationStats<UnderlyingStatsT>;
    std::vector<OptimizationStatsT> optimizationStatsVec;
    for (const SimplexTableauType simplexTableauType :
         absl::GetFlag(FLAGS_simplex_tableau_types)) {
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
          const std::string modelName = lpModelFileEntry.path().filename();
          SPDLOG_INFO("MODEL {}", modelName);
          auto linearProgram =
              MpsReader<FloatingPointT>().read(lpModelFileEntry.path());
          ASSERT_TRUE(linearProgram.has_value());

          if (isInstanceSuitable(basisSizeLimit, *linearProgram, modelName,
                                 isRelaxationOptType, allBoundsMustBeSpecified,
                                 instanceSetStats)) {
            lpOptFunc(*linearProgram, simplexTableauType,
                      lpModelFileEntry.path(), optimizationStatsVec);
          }
        }
      }
      SPDLOG_INFO("PROCESSED {} MODELS FROM {} TOTAL",
                  instanceSetStats._processedModelsCount,
                  instanceSetStats._totalModelsCount);
      SPDLOG_INFO("{} TOO BIG", instanceSetStats._tooBigModelsCount);
      SPDLOG_INFO("{} WITHOUT BOUNDS",
                  instanceSetStats._modelsWithNoBoundsCount);
      if (isRelaxationOptType) {
        ASSERT_EQ(instanceSetStats._totalModelsCount,
                  instanceSetStats._processedModelsCount +
                      instanceSetStats._tooBigModelsCount +
                      instanceSetStats._modelsWithNoBoundsCount);
      } else {
        SPDLOG_INFO("{} NOT PURE IP",
                    instanceSetStats._modelsWithNotAllIntegerVariables);
        SPDLOG_INFO("{} NOT ALL COEFFS INTEGER",
                    instanceSetStats._modelsWithNoAllIntegerCoeffs);
        SPDLOG_INFO("{} NOT ALL VARIABLES ARE NONNEGATIVE",
                    instanceSetStats._modelsWithNoAllNonnegativeVariables);
        ASSERT_EQ(instanceSetStats._totalModelsCount,
                  instanceSetStats._processedModelsCount +
                      instanceSetStats._tooBigModelsCount +
                      instanceSetStats._modelsWithNoBoundsCount +
                      instanceSetStats._modelsWithNotAllIntegerVariables +
                      instanceSetStats._modelsWithNoAllIntegerCoeffs +
                      instanceSetStats._modelsWithNoAllNonnegativeVariables);
      }
    }
    const bool extendedStatistics = absl::GetFlag(FLAGS_extended_statistics);
    if (!optimizationStatsVec.empty()) {
      OptStatisticsPrinter optStatisticsPrinter;
      optStatisticsPrinter.printFirstLine(
          optimizationStatsVec.front()._optimizationStats);
      for (const auto &optStats : optimizationStatsVec) {
        optStatisticsPrinter.print(optStats._optimizationStats,
                                   extendedStatistics,
                                   absl::GetFlag(FLAGS_cut_round_limit) / 10);
        optStatisticsPrinter.print(optStats._gurobiStats, false, 1);
      }
      SPDLOG_INFO(optStatisticsPrinter.toString());
    }
  }

  void compare(const LPOptimizationType lpOptimizationType,
               const LexicographicReoptType lexicographicReoptType,
               const IPOptStatistics<FloatingPointT> &ipOptStatistics,
               const LPOptStatistics<FloatingPointT> &gurobiLPOptStats,
               const isDualProgramOptimized isDualProgramOptimized =
                   isDualProgramOptimized::NO) {
    if (lpOptimizationType == LPOptimizationType::LINEAR_RELAXATION) {
      compareLPRelaxation(lexicographicReoptType, ipOptStatistics,
                          gurobiLPOptStats, isDualProgramOptimized);
    } else {
      compareIP(lexicographicReoptType, ipOptStatistics, gurobiLPOptStats);
    }
  }

  void
  compareLPRelaxation(const LexicographicReoptType lexicographicReoptType,
                      const IPOptStatistics<FloatingPointT> &ipOptStatistics,
                      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats,
                      const isDualProgramOptimized isDualProgramOptimized =
                          isDualProgramOptimized::NO) {
    const auto &firstRelaxationOptStats =
        ipOptStatistics._lpRelaxationStats.front();
    compareWithGurobi(lexicographicReoptType, firstRelaxationOptStats,
                      gurobiLPOptStats, isDualProgramOptimized);
  }
  void compareIP(const LexicographicReoptType lexicographicReoptType,
                 const IPOptStatistics<FloatingPointT> &ipOptStatistics,
                 const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    const auto &lastRelaxationOptStats =
        ipOptStatistics._lpRelaxationStats.back();
    compareWithGurobi(lastRelaxationOptStats._relaxationOptStats,
                      gurobiLPOptStats);
    checkIfSolutionIsInteger(ipOptStatistics);
  }

  void checkIfSolutionIsInteger(
      const IPOptStatistics<FloatingPointT> &ipOptStatistics) {
    const auto &optSolution = ipOptStatistics._optimalSolution;
    for (int solutionIdx = 0; solutionIdx < optSolution.size(); ++solutionIdx) {
      ASSERT_TRUE(NumericalTraitsT ::isInteger(optSolution[solutionIdx]))
          << fmt::format("Value {} at index {} is not integer",
                         optSolution[solutionIdx], solutionIdx);
    }
  }

  void compareWithGurobi(
      const LexicographicReoptType lexicographicReoptType,
      const LPRelaxationStatistics<FloatingPointT> &lpRelaxationStatistics,
      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats,
      const isDualProgramOptimized isDualProgramOptimized =
          isDualProgramOptimized::NO) {
    const auto &relaxationOptStats = lpRelaxationStatistics._relaxationOptStats;
    ASSERT_EQ(gurobiLPOptStats._optResult, relaxationOptStats._optResult);
    if (relaxationOptStats._optResult ==
        LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
      const auto optValueWithCorrectSign =
          isDualProgramOptimized == isDualProgramOptimized::YES
              ? -relaxationOptStats._optimalValue
              : relaxationOptStats._optimalValue;
      SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
      SPDLOG_INFO("SIMPLEX OPT AFTER INITIAL OPT {}", optValueWithCorrectSign);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue, optValueWithCorrectSign,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);

      const auto &optimalValueAfterFirstLexReopt =
          lpRelaxationStatistics._lexicographicReoptStats._optimalValue;
      const auto optValueAfterFirstLexReoptWithCorrectSign =
          isDualProgramOptimized == isDualProgramOptimized::YES
              ? -optimalValueAfterFirstLexReopt
              : optimalValueAfterFirstLexReopt;
      SPDLOG_INFO("SIMPLEX OPT AFTER LEXICOGRAPHIC {} REOPT {}",
                  lexicographicReoptTypeToStr(lexicographicReoptType),
                  optValueAfterFirstLexReoptWithCorrectSign);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  optValueAfterFirstLexReoptWithCorrectSign,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }

  void
  compareWithGurobi(const LPOptStatistics<FloatingPointT> &lpOptStatistics,
                    const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    ASSERT_EQ(gurobiLPOptStats._optResult, lpOptStatistics._optResult);
    if (lpOptStatistics._optResult ==
        LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
      SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
      SPDLOG_INFO("SIMPLEX OPT {}", lpOptStatistics._optimalValue);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue, lpOptStatistics._optimalValue,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }

  void compareWithGurobiDual(
      const LPOptStatistics<FloatingPointT> &primalProgramLpOptStatistics,
      const LPOptStatistics<FloatingPointT> &dualProgramLpOptStatistics,
      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    ASSERT_EQ(gurobiLPOptStats._optResult,
              primalProgramLpOptStatistics._optResult);
    ASSERT_EQ(gurobiLPOptStats._optResult,
              dualProgramLpOptStatistics._optResult);
    if (primalProgramLpOptStatistics._optResult ==
        LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
      SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
      SPDLOG_INFO("PRIMAL PROGRAM SIMPLEX OPT {}",
                  primalProgramLpOptStatistics._optimalValue);
      SPDLOG_INFO("DUAL PROGRAM SIMPLEX OPT {}",
                  -dualProgramLpOptStatistics._optimalValue);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  primalProgramLpOptStatistics._optimalValue,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  -dualProgramLpOptStatistics._optimalValue,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }

  void compareWithGurobiDual(
      const LPOptStatistics<FloatingPointT> &dualProgramLpOptStatistics,
      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    ASSERT_EQ(gurobiLPOptStats._optResult,
              dualProgramLpOptStatistics._optResult);
    if (dualProgramLpOptStatistics._optResult ==
        LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
      SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
      SPDLOG_INFO("DUAL PROGRAM SIMPLEX OPT {}",
                  -dualProgramLpOptStatistics._optimalValue);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  -dualProgramLpOptStatistics._optimalValue,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }

  void compareWithGurobi(
      const LexicographicReoptType lexicographicReoptType,
      const LexReoptStatistics<FloatingPointT> &lexReoptStatistics,
      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats,
      const std::vector<FloatingPointT> &gurobiSolution) {
    ASSERT_EQ(gurobiLPOptStats._optResult, lexReoptStatistics._optResult);
    if (lexReoptStatistics._optResult ==
        LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
      const auto &lexOptimizerSolution = lexReoptStatistics._solution;
      ASSERT_EQ(gurobiSolution.size(), lexOptimizerSolution.size());
      for (int varIdx = 0; varIdx < gurobiSolution.size(); ++varIdx) {
        EXPECT_NEAR(gurobiSolution[varIdx], lexOptimizerSolution[varIdx],
                    NumericalTraitsT::OPTIMALITY_TOLERANCE);
      }
      SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
      SPDLOG_INFO("LEXICOGRAPHIC {} SIMPLEX REOPT {}",
                  lexicographicReoptTypeToStr(lexicographicReoptType),
                  lexReoptStatistics._optimalValue);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  lexReoptStatistics._optimalValue,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }
};

#endif // GMISOLVER_LPTESTBASE_H
