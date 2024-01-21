#ifndef GMISOLVER_LEXREOPTSTATISTICS_H
#define GMISOLVER_LEXREOPTSTATISTICS_H

#include "src/DataModel/EnumTypes.h"
#include "src/Util/LPOptStatistics.h"

#include <numeric>
#include <string>
#include <vector>

template <typename T> struct LexReoptStatistics {
  using Type = T;

  void addLpOptStats(LPOptStatistics<T> lpOptStatistics) {
    _lexLPReoptStatsVec.push_back(std::move(lpOptStatistics));
  }

  size_t singleVarOptimizationCount() const {
    return _lexLPReoptStatsVec.size();
  }

  size_t totalIterationCount() const {
    return std::accumulate(
        _lexLPReoptStatsVec.begin(), _lexLPReoptStatsVec.end(), 0,
        [](const size_t value, const LPOptStatistics<T> &lpOptStatistics) {
          return value + lpOptStatistics._iterationCount;
        });
  }

  std::string _lpName;
  std::string _algorithmType;
  std::vector<LPOptStatistics<T>> _lexLPReoptStatsVec;
  std::vector<T> _solution;
  T _optimalValue;
  LPOptimizationResult _optResult;
  double _elapsedTimeSec{0.0};
  int32_t _reinversionFrequency;
};

#endif // GMISOLVER_LEXREOPTSTATISTICS_H
