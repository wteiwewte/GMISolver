#ifndef GMISOLVER_IPOPTSTATISTICS_H
#define GMISOLVER_IPOPTSTATISTICS_H

#include "src/Util/LPOptStatistics.h"

#include <map>
#include <string>
#include <vector>

template <typename T> struct LPRelaxationStatistics {
  LPOptStatistics<T> _relaxationOptStats;
  LexReoptStatistics<T> _lexicographicReoptStats;
};

template <typename T> struct IPOptStatistics {
  std::string _lpName;
  std::string _algorithmType;
  T _optimalValue{};
  std::vector<T> _optimalSolution;
  double _elapsedTimeMs{0.0};
  std::vector<LPRelaxationStatistics<T>> _lpRelaxationStats;
  int32_t _reinversionFrequency;
};

template <typename T>
using IPStatisticsMap = std::map<std::string, IPOptStatistics<T>>;

#endif // GMISOLVER_IPOPTSTATISTICS_H
