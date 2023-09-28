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
using LPOptStatisticsVec = std::vector<LPOptStatistics<T>>;

template <typename T> struct LexReoptStatistics {
  std::vector<LPOptStatistics<T>> _lexLPReoptStatsVec;
  std::vector<T> _solution;
  T _objectiveValueAfterLexReopt;
  LPOptimizationResult _optResult;
};

#endif // GMISOLVER_LPOPTSTATISTICS_H
