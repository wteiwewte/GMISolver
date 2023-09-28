#ifndef GMISOLVER_LPTESTBASE_H
#define GMISOLVER_LPTESTBASE_H

#include "src/Algorithms/SimplexTableau.h"
#include "src/DataModel/EnumTypes.h"
#include "src/Util/IPOptStatistics.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SpdlogHeader.h"

#include <tuple>

#include <gtest/gtest.h>

struct InstanceSetStats {
  int _processedModelsCount = 0;
  int _tooBigModelsCount = 0;
  int _modelsWithNoBoundsCount = 0;
  int _modelsWithNoAllIntegerVariables = 0;
  int _modelsWithNoAllIntegerCoeffs = 0;
  int _modelsWithNoAllNonnegativeVariables = 0;
  int _totalModelsCount = 0;
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
    //    if (modelName == "noswot.mps")
    //    {
    //      return false;
    //    }
    ++instanceSetStats._totalModelsCount;

    if (linearProgram.getRowInfos().size() > basisSizeLimit ||
        linearProgram.getVariableInfos().size() > 2 * basisSizeLimit) {
      ++instanceSetStats._tooBigModelsCount;
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

  template <typename LPOptFunc>
  void solveAndCompareInstancesFromSets(
      const std::string &dirPath, const size_t basisSizeLimit,
      const LPOptimizationType lpOptimizationType, LPOptFunc lpOptFunc,
      const bool allBoundsMustBeSpecified = true) {
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
            MpsReader<FloatingPointT>::read(lpModelFileEntry.path());
        ASSERT_TRUE(linearProgram.has_value());

        if (isInstanceSuitable(basisSizeLimit, *linearProgram, modelName,
                               isRelaxationOptType, allBoundsMustBeSpecified,
                               instanceSetStats)) {
          lpOptFunc(*linearProgram, lpModelFileEntry.path());
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
      SPDLOG_INFO("SIMPLEX OPT AFTER LEXICOGRAPHIC {} REOPT {}",
                  lexicographicReoptTypeToStr(lexicographicReoptType),
                  optimalValueAfterFirstLexReopt);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  optimalValueAfterFirstLexReopt,
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
      SPDLOG_INFO("SIMPLEX OPT AFTER OPT {}", lpOptStatistics._optimalValue);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue, lpOptStatistics._optimalValue,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }
  void compareWithGurobi(
      const LexicographicReoptType lexicographicReoptType,
      const LexReoptStatistics<FloatingPointT> &lexReoptStatistics,
      const LPOptStatistics<FloatingPointT> &gurobiLPOptStats) {
    ASSERT_EQ(gurobiLPOptStats._optResult, lexReoptStatistics._optResult);
    if (lexReoptStatistics._optResult ==
        LPOptimizationResult::BOUNDED_AND_FEASIBLE) {
      SPDLOG_INFO("GUROBI OPT {}", gurobiLPOptStats._optimalValue);
      SPDLOG_INFO("SIMPLEX OPT AFTER LEXICOGRAPHIC {} REOPT {}",
                  lexicographicReoptTypeToStr(lexicographicReoptType),
                  lexReoptStatistics._objectiveValueAfterLexReopt);
      EXPECT_NEAR(gurobiLPOptStats._optimalValue,
                  lexReoptStatistics._objectiveValueAfterLexReopt,
                  NumericalTraitsT::OPTIMALITY_TOLERANCE);
    }
  }
};

#endif // GMISOLVER_LPTESTBASE_H
