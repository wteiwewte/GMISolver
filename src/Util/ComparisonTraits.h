#ifndef GMISOLVER_COMPARISONTRAITS_H
#define GMISOLVER_COMPARISONTRAITS_H

#include <algorithm>
#include <cmath>
#include <numeric>

template <typename T> struct ApproximateComparisonTraits {
  static_assert(!std::numeric_limits<T>::is_integer);

  constexpr static T ABSOLUTE_EPSILON = 1e-9;
  constexpr static T PIVOT_ABSOLUTE_EPSILON = 1e-6;
  //  constexpr static T ABSOLUTE_EPSILON = std::numeric_limits<T>::epsilon();
  constexpr static T RELATIVE_EPSILON = std::numeric_limits<T>::epsilon();
  constexpr static size_t ULP_COUNT = 4;

  static bool equal(const T &x, const T &y) {
    const T diff = std::fabs(x - y);
    return (diff < ABSOLUTE_EPSILON) ||
           (diff <= RELATIVE_EPSILON * std::fabs(x + y) * ULP_COUNT);
  }

  static bool greater(const T &x, const T &y) {
    return (x > y) && !equal(x, y);
  }

  static bool less(const T &x, const T &y) { return (x < y) && !equal(x, y); }

  static bool isEligibleForPivot(const T &x) {
    return std::fabs(x) > PIVOT_ABSOLUTE_EPSILON;
  }
  //  static bool canBeOmittedInComputations(const T &x) { return std::fabs(x) <
  //  PIVOT_ABSOLUTE_EPSILON; }
};

template <typename T> struct SimpleComparisonTraits {
  static bool equal(const T &x, const T &y) { return x == y; }

  static bool less(const T &x, const T &y) { return x < y; }

  static bool greater(const T &x, const T &y) { return x > y; }
};

#endif // GMISOLVER_COMPARISONTRAITS_H
