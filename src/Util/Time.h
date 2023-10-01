#ifndef GMISOLVER_TIME_H
#define GMISOLVER_TIME_H

#include <chrono>

template <typename Func, typename TimeRep = std::chrono::microseconds>
double executeAndMeasureTime(Func f) {
  auto start = std::chrono::steady_clock::now();
  f();
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end - start)
             .count() /
         1000000.0;
}

#endif // GMISOLVER_TIME_H
