#ifndef GMISOLVER_LPOPTSTATISTICS_H
#define GMISOLVER_LPOPTSTATISTICS_H

#include "src/DataModel/EnumTypes.h"

#include <string>
#include <vector>

template <typename T> struct LPOptStatistics {
  using Type = T;

  std::string _lpName;
  std::string _algorithmType;
  std::vector<T> _consecutiveObjectiveValues;
  LPOptimizationResult _optResult;
  T _optimalValue{};
  int32_t _iterationCount{0};
  double _elapsedTimeSec{0.0};
  bool _phaseOneSucceeded{false}; // only for primal simplex
  int32_t _reinversionFrequency;
};

template <typename T>
using LPOptStatisticsVec = std::vector<LPOptStatistics<T>>;

template <typename T> struct PrimalSimplexOutput {
  using Type = T;

  LPOptStatistics<T> _phaseOneLpOptStats;
  std::optional<LPOptStatistics<T>> _phaseTwoLpOptStats;
};

#endif // GMISOLVER_LPOPTSTATISTICS_H
