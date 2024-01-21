#ifndef GMISOLVER_OPTIMIZATIONSTATSPRINTER_H
#define GMISOLVER_OPTIMIZATIONSTATSPRINTER_H

#include "src/Util/IPOptStatistics.h"
#include "src/Util/LPOptStatistics.h"
#include "src/Util/LexReoptStatistics.h"

#include <iomanip>
#include <map>
// fmt must be included after spdlog - weird bug
#include <fmt/format.h>
#include <sstream>
#include <string>

struct OptStatisticsPrinter {
  constexpr static int LP_NAME_WIDTH = 30;
  constexpr static int ALGO_TYPE_WIDTH = 71;
  constexpr static int OPT_RESULT_WIDTH = 25;
  constexpr static int OPT_VALUE_WIDTH = 20;
  constexpr static int OPT_COUNT_WIDTH = 11;
  constexpr static int ITERATION_COUNT_WIDTH = 12;
  constexpr static int ELAPSED_TIME_WIDTH = 20;
  constexpr static int REINVERSION_FREQUENCY_WIDTH = 12;
  constexpr static int FLOAT_NUMBERS_PRECISION = 6;

  OptStatisticsPrinter() {
    _totalWidth = LP_NAME_WIDTH + ALGO_TYPE_WIDTH + OPT_RESULT_WIDTH +
                  OPT_VALUE_WIDTH + ITERATION_COUNT_WIDTH +
                  REINVERSION_FREQUENCY_WIDTH;

    _oss << '\n';
  }

  void printFirstLineBasicParameters() {
    _oss << fmt::format("{:^{}}|", "LP NAME", LP_NAME_WIDTH);
    _oss << fmt::format("{:^{}}|", "ALGO TYPE", ALGO_TYPE_WIDTH);
    _oss << fmt::format("{:^{}}|", "OPT RESULT", OPT_RESULT_WIDTH);
    _oss << fmt::format("{:^{}}|", "OPT VALUE", OPT_VALUE_WIDTH);
    _oss << fmt::format("{:^{}}|", "ITER COUNT", ITERATION_COUNT_WIDTH);
    _oss << fmt::format("{:^{}}|", "ELP TIME SECONDS", ELAPSED_TIME_WIDTH);
    _oss << fmt::format("{:^{}}|", "REINV FREQ", REINVERSION_FREQUENCY_WIDTH);
  }

  template <typename T> void printFirstLine(const LPOptStatistics<T> &) {
    printFirstLineBasicParameters();
    printLineBreak();
  }

  template <typename T>
  void print(const LPOptStatistics<T> &lpOptStatistics, const bool) {
    _oss << fmt::format("{:^{}}|", lpOptStatistics._lpName, LP_NAME_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptStatistics._algorithmType,
                        ALGO_TYPE_WIDTH);
    _oss << fmt::format("{:^{}}|",
                        lpOptimizationResultToStr(lpOptStatistics._optResult),
                        OPT_RESULT_WIDTH);
    _oss << valueWithPrecision(lpOptStatistics._optimalValue, OPT_VALUE_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptStatistics._iterationCount,
                        ITERATION_COUNT_WIDTH);
    _oss << valueWithPrecision(lpOptStatistics._elapsedTimeSec,
                               ELAPSED_TIME_WIDTH);
    _oss << fmt::format("{:^{}}|", lpOptStatistics._reinversionFrequency,
                        REINVERSION_FREQUENCY_WIDTH);
    printLineBreak();
  }

  template <typename T> void printFirstLine(const IPOptStatistics<T> &) {
    printFirstLineBasicParameters();
    _oss << fmt::format("{:^{}}|", "RELAX OPT COUNT", OPT_COUNT_WIDTH);
    printLineBreak();
  }

  template <typename T>
  void print(const IPOptStatistics<T> &ipOptStatistics, const bool extended) {
    _oss << fmt::format("{:^{}}|", ipOptStatistics._lpName, LP_NAME_WIDTH);
    _oss << fmt::format("{:^{}}|", ipOptStatistics._algorithmType,
                        ALGO_TYPE_WIDTH);
    _oss << fmt::format("{:^{}}|",
                        lpOptimizationResultToStr(ipOptStatistics._optResult),
                        OPT_RESULT_WIDTH);
    _oss << valueWithPrecision(ipOptStatistics._optimalValue, OPT_VALUE_WIDTH);
    _oss << fmt::format("{:^{}}|", ipOptStatistics.totalIterationCount(),
                        ITERATION_COUNT_WIDTH);
    _oss << valueWithPrecision(ipOptStatistics._elapsedTimeSec,
                               ELAPSED_TIME_WIDTH);
    _oss << fmt::format("{:^{}}|", ipOptStatistics._reinversionFrequency,
                        REINVERSION_FREQUENCY_WIDTH);
    _oss << fmt::format("{:^{}}|",
                        ipOptStatistics.relaxationOptimizationCount(),
                        OPT_COUNT_WIDTH);
    printLineBreak();

    if (extended) {
      for (const auto &lpRelaxationStats : ipOptStatistics._lpRelaxationStats) {
        print(lpRelaxationStats._relaxationOptStats, false);
        print(lpRelaxationStats._lexicographicReoptStats, false);
      }
    }
  }

  template <typename T> void printFirstLine(const LexReoptStatistics<T> &) {
    printFirstLineBasicParameters();
    _oss << fmt::format("{:^{}}|", "LP OPT COUNT", OPT_COUNT_WIDTH);
    printLineBreak();
  }

  template <typename T>
  void print(const LexReoptStatistics<T> &lexReoptStatistics,
             const bool extended) {
    _oss << fmt::format("{:^{}}|", lexReoptStatistics._lpName, LP_NAME_WIDTH);
    _oss << fmt::format("{:^{}}|", lexReoptStatistics._algorithmType,
                        ALGO_TYPE_WIDTH);
    _oss << fmt::format(
        "{:^{}}|", lpOptimizationResultToStr(lexReoptStatistics._optResult),
        OPT_RESULT_WIDTH);
    _oss << valueWithPrecision(lexReoptStatistics._optimalValue,
                               OPT_VALUE_WIDTH);
    _oss << fmt::format("{:^{}}|", lexReoptStatistics.totalIterationCount(),
                        ITERATION_COUNT_WIDTH);
    _oss << valueWithPrecision(lexReoptStatistics._elapsedTimeSec,
                               ELAPSED_TIME_WIDTH);
    _oss << fmt::format("{:^{}}|", lexReoptStatistics._reinversionFrequency,
                        REINVERSION_FREQUENCY_WIDTH);
    _oss << fmt::format("{:^{}}|",
                        lexReoptStatistics.singleVarOptimizationCount(),
                        OPT_COUNT_WIDTH);
    printLineBreak();

    if (extended) {
      for (const auto &lpOpt : lexReoptStatistics._lexLPReoptStatsVec) {
        print(lpOpt, false);
      }
    }
  }

  template <typename T>
  void printFirstLine(const PrimalSimplexOutput<T> &primalSimplexOutput) {
    printFirstLine(primalSimplexOutput._phaseOneLpOptStats);
  }

  template <typename T>
  void print(const PrimalSimplexOutput<T> &primalSimplexOutput,
             const bool extended) {
    print(primalSimplexOutput._phaseOneLpOptStats, extended);
    if (primalSimplexOutput._phaseTwoLpOptStats.has_value()) {
      print(primalSimplexOutput._phaseTwoLpOptStats.value(), extended);
    }
  }

  template <typename T>
  std::string valueWithPrecision(const T value, const int precision) const {
    const auto formattedValue = fmt::format(
        "{:^{}}", fmt::format("{:.{}f}", value, FLOAT_NUMBERS_PRECISION),
        precision);
    const auto dotPosition = formattedValue.find('.');
    constexpr static int FIXED_DOT_POSITION = 7;

    auto finalResult = formattedValue;
    if (dotPosition < FIXED_DOT_POSITION) {
      finalResult =
          std::string(FIXED_DOT_POSITION - dotPosition, ' ') + formattedValue;
      finalResult = finalResult.substr(0, formattedValue.size());
    }

    if (dotPosition > FIXED_DOT_POSITION) {
      finalResult = formattedValue.substr(dotPosition - FIXED_DOT_POSITION);
      finalResult += std::string(dotPosition - FIXED_DOT_POSITION, ' ');
    }

    return finalResult + '|';
  }

  void printLineBreak() {
    _oss.width(_totalWidth);
    _oss.fill('-');
    _oss << '-' << '\n';
    _oss.fill(' ');
  }

  std::string toString() const { return _oss.str(); }

  std::ostringstream _oss;
  int _totalWidth{0};
};

#endif // GMISOLVER_OPTIMIZATIONSTATSPRINTER_H
