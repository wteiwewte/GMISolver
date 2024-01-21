#ifndef GMISOLVER_IPOPTSTATISTICS_H
#define GMISOLVER_IPOPTSTATISTICS_H

#include "src/Util/LPOptStatistics.h"
#include "src/Util/LexReoptStatistics.h"

#include <map>
#include <string>
#include <vector>

template <typename T> struct LPRelaxationStatistics {
  LPOptStatistics<T> _relaxationOptStats;
  LexReoptStatistics<T> _lexicographicReoptStats;
};

template <typename T> struct IPOptStatistics {
  using Type = T;

  size_t relaxationOptimizationCount() const {
    return _lpRelaxationStats.size();
  }

  size_t totalIterationCount() const {
    return std::accumulate(
        _lpRelaxationStats.begin(), _lpRelaxationStats.end(), 0,
        [](const size_t value,
           const LPRelaxationStatistics<T> &lpRelaxationStatistics) {
          return value +
                 lpRelaxationStatistics._relaxationOptStats._iterationCount +
                 lpRelaxationStatistics._lexicographicReoptStats
                     .totalIterationCount();
        });
  }

  std::string _lpName;
  std::string _algorithmType;
  T _optimalValue{};
  std::vector<T> _optimalSolution;
  LPOptimizationResult _optResult;
  double _elapsedTimeSec{0.0};
  std::vector<LPRelaxationStatistics<T>> _lpRelaxationStats;
  int32_t _reinversionFrequency;
};

template <typename T>
using IPStatisticsMap = std::map<std::string, IPOptStatistics<T>>;

#endif // GMISOLVER_IPOPTSTATISTICS_H
