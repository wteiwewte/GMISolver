#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/LPOptStatisticsPrinter.h"
#include "src/Util/MpsReader.h"
#include "src/Util/SpdlogHeader.h"

#include <filesystem>
#include <optional>
#include <thread>
#include <vector>

#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

ABSL_FLAG(std::string, log_level, "info", "log level");
ABSL_FLAG(std::string, log_file, "GMISolver_out.txt", "log file");
ABSL_FLAG(std::string, lp_opt_stats_file, "LPOptStatistics.txt", "log file");
ABSL_FLAG(std::vector<std::string>, benchmark_simplex_types,
          std::vector<std::string>({"primal", "dual"}),
          "List of simplex algorithms to benchmark");
ABSL_FLAG(std::optional<std::string>, lp_model_file, std::nullopt,
          "Path to single lp model");
ABSL_FLAG(std::optional<std::string>, lp_models_directory, std::nullopt,
          "Path to directory with lp models");
ABSL_FLAG(
    int32_t, obj_value_logging_frequency, 100,
    "Current objective value should be logged every nth iteration of simplex");
ABSL_FLAG(int32_t, reinversion_frequency, 300,
          "Basis matrix should be reinverted every nth iteration of simplex");
ABSL_FLAG(bool, use_product_form_of_inverse, true,
          "Basis matrix inverse is represented via product form of inverse");

bool contains(const std::vector<std::string> &vec, const std::string &str) {
  return std::find(vec.begin(), vec.end(), str) != vec.end();
}

template <typename T, typename SimplexTraitsT>
void runPrimalSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram, LPOptStatisticsVec<T>& lpOptStatisticsVec) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, true, absl::GetFlag(FLAGS_use_product_form_of_inverse));
  RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT> revisedPrimalSimplexPfiBounds(
      simplexTableau,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_reinversion_frequency));
  auto phaseOneLpOptStats = revisedPrimalSimplexPfiBounds.runPhaseOne();
  lpOptStatisticsVec.push_back(phaseOneLpOptStats);
  if (!phaseOneLpOptStats._phaseOneSucceeded) {
    SPDLOG_WARN("PHASE ONE OF {} ALGORITHM FAILED", revisedPrimalSimplexPfiBounds.name());
    return;
  }
  lpOptStatisticsVec.push_back(revisedPrimalSimplexPfiBounds.runPhaseTwo());

  SPDLOG_INFO(simplexTableau.toStringObjectiveValue());
//  SPDLOG_INFO(simplexTableau.toStringSolution());
//  revisedPrimalSimplexPfiBounds.lexicographicReoptimization(true);
//  SPDLOG_INFO(simplexTableau.toStringObjectiveValue());
//  SPDLOG_INFO(simplexTableau.toStringSolution());
}

template <typename T, typename SimplexTraitsT>
void runDualSimplexWithImplicitBounds(const LinearProgram<T> &linearProgram, LPOptStatisticsVec<T>& lpOptStatisticsVec) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, false, absl::GetFlag(FLAGS_use_product_form_of_inverse));
  lpOptStatisticsVec.push_back(RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
      simplexTableau, DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_reinversion_frequency))
      .run());
  SPDLOG_INFO(simplexTableau.toStringObjectiveValue());
//  SPDLOG_INFO(simplexTableau.toStringSolution());

  // primal lex opt
//  RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT> revisedPrimalSimplexPfiBounds(
//      simplexTableau,
//      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
//      absl::GetFlag(FLAGS_obj_value_logging_frequency),
//      absl::GetFlag(FLAGS_reinversion_frequency));
//  revisedPrimalSimplexPfiBounds.lexicographicReoptimization(false);
//  SPDLOG_INFO(simplexTableau.toStringObjectiveValue());
//  SPDLOG_INFO(simplexTableau.toStringSolution());
}

template <typename T>
bool isLPTooBig(const LinearProgram<T>& linearProgram)
{
  constexpr int MAX_NUMBER_OF_CONTRAINTS = 500;
  constexpr int MAX_NUMBER_OF_VARIABLES = 500;
  return linearProgram.getRowInfos().size() > MAX_NUMBER_OF_CONTRAINTS || linearProgram.getVariableInfos().size() > MAX_NUMBER_OF_VARIABLES;
}

template <typename T>
void readLPModelAndProcess(const std::filesystem::path &modelFileMpsPath, LPOptStatisticsVec<T>& lpOptStatisticsVec) {
  SPDLOG_INFO("Processing lp model from file {}", std::string{modelFileMpsPath});
  auto linearProgram = MpsReader<T>::read(modelFileMpsPath);
  if (!linearProgram.has_value()) {
    spdlog::warn("Could not read properly lp from {} file", std::string{modelFileMpsPath});
    lpOptStatisticsVec.push_back(LPOptStatistics<T>{._lpName = modelFileMpsPath.filename()});
    return;
  }

  constexpr bool SKIP_TO_BIG_LPS = false;
  if (SKIP_TO_BIG_LPS && isLPTooBig(*linearProgram))
  {
    SPDLOG_INFO("Skipped too big LP from file {}", std::string{modelFileMpsPath});
    return;
  }

  const auto simplexTypesToBenchmark =
      absl::GetFlag(FLAGS_benchmark_simplex_types);
  const bool processAllTypes = contains(simplexTypesToBenchmark, "all");
  if (processAllTypes || contains(simplexTypesToBenchmark, "primal"))
    runPrimalSimplexWithImplicitBounds<T, SimplexTraits<T, false>>(*linearProgram, lpOptStatisticsVec);

  if (processAllTypes || contains(simplexTypesToBenchmark, "primal_sparse"))
    runPrimalSimplexWithImplicitBounds<T, SimplexTraits<T, true>>(*linearProgram, lpOptStatisticsVec);

  if (processAllTypes || contains(simplexTypesToBenchmark, "dual"))
    runDualSimplexWithImplicitBounds<T, SimplexTraits<T, false>>(*linearProgram, lpOptStatisticsVec);

  if (processAllTypes || contains(simplexTypesToBenchmark, "dual_sparse"))
    runDualSimplexWithImplicitBounds<T, SimplexTraits<T, true>>(*linearProgram, lpOptStatisticsVec);
}

void initFileLogger() {
  auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_st>(
      absl::GetFlag(FLAGS_log_file), true);
  auto fileLogger = std::make_shared<spdlog::logger>("fileLogger", fileSink);
  spdlog::set_default_logger(fileLogger);

  const std::string COUT_PATTERN = "%v";
  spdlog::set_pattern("%^[%Y-%m-%d %H:%M:%S.%e][%L]%$ %v");
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_every(std::chrono::seconds(5));
}

template <typename T>
void printLPOptStats(const LPOptStatisticsVec<T>& lpOptStatisticsVec)
{
  auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_st>(
      absl::GetFlag(FLAGS_lp_opt_stats_file), true);
  auto fileLogger = std::make_shared<spdlog::logger>("fileLogger", fileSink);
  LPOptStatisticsPrinter<T> lpOptStatisticsPrinter(lpOptStatisticsVec);
  lpOptStatisticsPrinter.print();
  SPDLOG_LOGGER_INFO(fileLogger, lpOptStatisticsPrinter.toString());
}

int main(int argc, char **argv) {
  absl::ParseCommandLine(argc, argv);
  initFileLogger();
  std::this_thread::sleep_for(std::chrono::seconds(1));

  using FloatingPointT = double;

  LPOptStatisticsVec<FloatingPointT> lpOptStatisticsVec;

  if (const auto lpModelFile = absl::GetFlag(FLAGS_lp_model_file);
      lpModelFile.has_value())
    readLPModelAndProcess<FloatingPointT>(std::filesystem::path{*lpModelFile}, lpOptStatisticsVec);
  else if (const auto lpModelsDirectory =
               absl::GetFlag(FLAGS_lp_models_directory);
           lpModelsDirectory.has_value())
    for (const auto &lpModelFileEntry :
         std::filesystem::directory_iterator(*lpModelsDirectory))
      readLPModelAndProcess<FloatingPointT>(lpModelFileEntry.path(), lpOptStatisticsVec);

  printLPOptStats(lpOptStatisticsVec);
}
