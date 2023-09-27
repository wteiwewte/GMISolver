#ifndef GMISOLVER_LPOPTSTATISTICS_H
#define GMISOLVER_LPOPTSTATISTICS_H

#include "src/DataModel/EnumTypes.h"

#include <string>
#include <vector>

template <typename T> struct LPOptStatistics {
  std::string _lpName;
  std::string _simplexAlgorithmType;
  std::vector<T> _consecutiveObjectiveValues;
  LPOptimizationResult _optResult;
  T _optimalValue{};
  int32_t _iterationCount{0};
  bool _phaseOneSucceeded{false}; // only for primal simplex
  int32_t _reinversionFrequency;
};

template <typename T>
LPOptStatistics<T>
convert(const LPOptStatistics<double> &doubleLpOptStatistics) {
  return LPOptStatistics<T>{
      ._lpName = doubleLpOptStatistics._lpName,
      ._simplexAlgorithmType = doubleLpOptStatistics._simplexAlgorithmType,
      ._consecutiveObjectiveValues =
          std::vector<T>{
              doubleLpOptStatistics._consecutiveObjectiveValues.begin(),
              doubleLpOptStatistics._consecutiveObjectiveValues.end()},
      ._optResult = doubleLpOptStatistics._optResult,
      ._optimalValue = T{doubleLpOptStatistics._optimalValue},
      ._iterationCount = doubleLpOptStatistics._iterationCount,
      ._phaseOneSucceeded = doubleLpOptStatistics._phaseOneSucceeded,
      ._reinversionFrequency = doubleLpOptStatistics._reinversionFrequency,
  };
}

template <typename T>
using LPOptStatisticsVec = std::vector<LPOptStatistics<T>>;

template <typename T> struct LexReoptStatistics {
  std::vector<LPOptStatistics<T>> _lexLPReoptStatsVec;
  T _objectiveValueAfterLexReopt;
  LPOptimizationResult _optResult;
};

#endif // GMISOLVER_LPOPTSTATISTICS_H
