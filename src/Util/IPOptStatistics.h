#ifndef GMISOLVER_LPOPTSTATISTICS_H
#define GMISOLVER_LPSTATISTICS_H

#include "src/Util/LPOptStatistics.h"

#include <map>
#include <string>
#include <vector>

template <typename T> struct LPRelaxationStatistics {
  LPOptStatistics<T> _relaxationOptStats;
  LexReoptStatistics<T> _lexicographicReoptStats;
};

template <typename T> struct IPOptStatistics {
  T _optimalValue{};
  std::vector<T> _optimalSolution;
  std::vector<LPRelaxationStatistics<T>> _lpRelaxationStats;
};

template <typename T>
using IPStatisticsMap = std::map<std::string, IPOptStatistics<T>>;

#endif // GMISOLVER_LPOPTSTATISTICS_H
