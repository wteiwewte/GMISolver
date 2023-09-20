#ifndef GMISOLVER_LPOPTSTATISTICS_H
#define GMISOLVER_LPSTATISTICS_H

#include "src/Util/LPOptStatistics.h"

#include <map>
#include <string>
#include <vector>

template <typename T> struct IPOptStatistics {
  T _integerOptimum{};
  std::vector<LPStatistics<T>> _lpRelaxationsStats;
};

template <typename T>
using IPStatisticsMap = std::map<std::string, IPOptStatistics<T>>;

#endif // GMISOLVER_LPOPTSTATISTICS_H
