#ifndef GMISOLVER_LPOPTSTATISTICSPRINTER_H
#define GMISOLVER_LPOPTSTATISTICSPRINTER_H

#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/MatrixTypes.h"
#include "src/DataModel/SimplexBasisData.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/SpdlogHeader.h"

#include <iomanip>
#include <map>
// fmt must be included after spdlog - weird bug
#include <fmt/format.h>
#include <functional>
#include <optional>
#include <sstream>
#include <string>

template <typename T>
struct LPOptStatisticsPrinter {
  constexpr static int LP_NAME_WIDTH = 20;
  constexpr static int SIMPLEX_ALGO_TYPE_WIDTH = 35;
  constexpr static int LP_OPT_RESULT_WIDTH = 25;
  constexpr static int OPT_VALUE_WIDTH = 20;
  constexpr static int ITERATION_COUNT_WIDTH = 12;
  constexpr static int REINVERSION_FREQUENCY_WIDTH = 12;


  LPOptStatisticsPrinter(const LPOptStatisticsVec<T>& lpOptStatisticsVec)
      : _lpOptStatisticsVec(lpOptStatisticsVec) {
    _totalWidth =
        LP_NAME_WIDTH + SIMPLEX_ALGO_TYPE_WIDTH + LP_OPT_RESULT_WIDTH + OPT_VALUE_WIDTH + ITERATION_COUNT_WIDTH + REINVERSION_FREQUENCY_WIDTH;

    _oss << '\n';
  }

  void print()
  {
    printFirstLine();
    for (const auto& lpOptStatistics : _lpOptStatisticsVec)
    {
      printLPOptStatistics(lpOptStatistics);
    }
  }

  void printFirstLine()
  {
    _oss << fmt::format("{:^{}}|", "LP NAME", LP_NAME_WIDTH);
    _oss << fmt::format("{:^{}}|", "SIMPLEX ALGO TYPE", SIMPLEX_ALGO_TYPE_WIDTH);
    _oss << fmt::format("{:^{}}|", "LP OPT RESULT", LP_OPT_RESULT_WIDTH);
    _oss << fmt::format("{:^{}}|", "OPT VALUE", OPT_VALUE_WIDTH);
    _oss << fmt::format("{:^{}}|", "ITER COUNT", ITERATION_COUNT_WIDTH);
    _oss << fmt::format("{:^{}}|", "REINV FREQ", REINVERSION_FREQUENCY_WIDTH);
    printLineBreak();
  }

  void printLPOptStatistics(const LPOptStatistics<T>& lpOptStatistics)
  {
    _oss << fmt::format("{:^{}}|", lpOptStatistics._lpName, LP_NAME_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptStatistics._simplexAlgorithmType, SIMPLEX_ALGO_TYPE_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptimizationResultToStr(lpOptStatistics._optResult), LP_OPT_RESULT_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptStatistics._optimalValue, OPT_VALUE_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptStatistics._iterationCount, ITERATION_COUNT_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptStatistics._reinversionFrequency, REINVERSION_FREQUENCY_WIDTH);
    printLineBreak();
  }

  void printLineBreak() {
    _oss.width(_totalWidth);
    _oss.fill('-');
    _oss << '-' << '\n';
    _oss.fill(' ');
  }


  std::string toString() const { return _oss.str(); }

  const LPOptStatisticsVec<T>& _lpOptStatisticsVec;
  std::ostringstream _oss;
  int _totalWidth{0};
};

#endif // GMISOLVER_LPOPTSTATISTICSPRINTER_H
