#include "src/Algorithms/DualSimplex.h"
#include "src/Algorithms/DualSimplexGomory.h"
#include "src/Algorithms/PrimalSimplex.h"
#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/GurobiOptimizer.h"
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
ABSL_FLAG(std::string, lp_opt_stats_file, "LPOptStatistics.txt",
          "Lp opt statistics log file");
ABSL_FLAG(std::string, gurobi_log_file, "Gurobi_logs.txt", "Gurobi log file");
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
ABSL_FLAG(std::vector<SimplexTableauType>, simplex_tableau_types,
          {SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE},
          "Simplex tableau types");
ABSL_FLAG(ValidateSimplexOption, validate_simplex_option,
          ValidateSimplexOption::DONT_VALIDATE,
          "Validate simplex implementations");
ABSL_FLAG(bool, run_gomory, false, "Run gomory scheme");

bool contains(const std::vector<std::string> &vec, const std::string &str) {
  return std::find(vec.begin(), vec.end(), str) != vec.end();
}

template <typename T, typename SimplexTraitsT>
void runPrimalSimplexWithImplicitBounds(
    const LinearProgram<T> &linearProgram,
    const SimplexTableauType simplexTableauType,
    LPOptStatisticsVec<T> &lpOptStatisticsVec) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::PRIMAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  PrimalSimplex<T, SimplexTraitsT> revisedPrimalSimplexPfiBounds(
      simplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));
  auto phaseOneLpOptStats = revisedPrimalSimplexPfiBounds.runPhaseOne();
  lpOptStatisticsVec.push_back(phaseOneLpOptStats);
  if (!phaseOneLpOptStats._phaseOneSucceeded) {
    SPDLOG_WARN("PHASE ONE OF {} ALGORITHM FAILED",
                revisedPrimalSimplexPfiBounds.type());
    return;
  }
  lpOptStatisticsVec.push_back(revisedPrimalSimplexPfiBounds.runPhaseTwo());

  SPDLOG_INFO(simplexTableau.toStringObjectiveValue());
}

template <typename T, typename SimplexTraitsT>
void runDualSimplexWithImplicitBounds(
    const LinearProgram<T> &linearProgram,
    const SimplexTableauType simplexTableauType,
    LPOptStatisticsVec<T> &lpOptStatisticsVec) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  lpOptStatisticsVec.push_back(
      DualSimplex<T, SimplexTraitsT>(
          simplexTableau, reinversionManager,
          DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
          absl::GetFlag(FLAGS_obj_value_logging_frequency),
          absl::GetFlag(FLAGS_validate_simplex_option))
          .run(""));
  SPDLOG_INFO(simplexTableau.toStringObjectiveValue());
}

template <typename T, typename SimplexTraitsT>
void runDualSimplexGomoryWithPrimalCuts(
    const LinearProgram<T> &linearProgram,
    const SimplexTableauType simplexTableauType,
    LPOptStatisticsVec<T> &lpOptStatisticsVec) {
  SimplexTableau<T, SimplexTraitsT> simplexTableau(
      linearProgram, SimplexType::DUAL, simplexTableauType);
  ReinversionManager<T, SimplexTraitsT> reinversionManager(
      simplexTableau, absl::GetFlag(FLAGS_reinversion_frequency));
  DualSimplexGomory<T, SimplexTraitsT> dualSimplexGomoryWithPrimalCuts(
      simplexTableau, reinversionManager,
      PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST,
      DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION,
      absl::GetFlag(FLAGS_obj_value_logging_frequency),
      absl::GetFlag(FLAGS_validate_simplex_option));
  const auto ipOptStats = dualSimplexGomoryWithPrimalCuts.run(
      LexicographicReoptType::MAX, LPOptimizationType::LINEAR_RELAXATION,
      GomoryCutChoosingRule::FIRST);
  for (const auto &lpRelaxStats : ipOptStats._lpRelaxationStats) {
    lpOptStatisticsVec.push_back(lpRelaxStats._relaxationOptStats);
    const auto &lexLPReoptStatsVec =
        lpRelaxStats._lexicographicReoptStats._lexLPReoptStatsVec;
    lpOptStatisticsVec.insert(lpOptStatisticsVec.end(),
                              lexLPReoptStatsVec.begin(),
                              lexLPReoptStatsVec.end());
  }
}

template <typename T> bool isLPTooBig(const LinearProgram<T> &linearProgram) {
  constexpr int MAX_NUMBER_OF_CONTRAINTS = 500;
  constexpr int MAX_NUMBER_OF_VARIABLES = 500;
  return linearProgram.getRowInfos().size() > MAX_NUMBER_OF_CONTRAINTS ||
         linearProgram.getVariableInfos().size() > MAX_NUMBER_OF_VARIABLES;
}

template <typename T>
void readLPModelAndProcess(const std::filesystem::path &modelFileMpsPath,
                           const SimplexTableauType simplexTableauType,
                           LPOptStatisticsVec<T> &lpOptStatisticsVec) {
  SPDLOG_INFO("Processing lp model from file {}",
              std::string{modelFileMpsPath});
  auto linearProgram = MpsReader<T>::read(modelFileMpsPath);
  if (!linearProgram.has_value()) {
    spdlog::warn("Could not read properly lp from {} file",
                 std::string{modelFileMpsPath});
    lpOptStatisticsVec.push_back(
        LPOptStatistics<T>{._lpName = modelFileMpsPath.filename(),
                           ._optResult = LPOptimizationResult::COULD_NOT_LOAD});
    return;
  }

  constexpr bool SKIP_TO_BIG_LPS = false;
  if (SKIP_TO_BIG_LPS && isLPTooBig(*linearProgram)) {
    SPDLOG_INFO("Skipped too big LP from file {}",
                std::string{modelFileMpsPath});
    return;
  }
  if (!linearProgram->checkIfAllBoundsAreSpeficied()) {
    SPDLOG_INFO(
        "SKIPPING MODEL {} BECAUSE NOT ALL VARIABLES HAVE BOUNDS SPECIFIED",
        std::string{modelFileMpsPath});
    return;
  }

  const auto simplexTypesToBenchmark =
      absl::GetFlag(FLAGS_benchmark_simplex_types);
  const bool processAllTypes = contains(simplexTypesToBenchmark, "all");
  if (processAllTypes || contains(simplexTypesToBenchmark, "primal"))
    runPrimalSimplexWithImplicitBounds<
        T, SimplexTraits<T, MatrixRepresentationType::NORMAL>>(
        *linearProgram, simplexTableauType, lpOptStatisticsVec);

  if (processAllTypes || contains(simplexTypesToBenchmark, "primal_sparse"))
    runPrimalSimplexWithImplicitBounds<
        T, SimplexTraits<T, MatrixRepresentationType::SPARSE>>(
        *linearProgram, simplexTableauType, lpOptStatisticsVec);

  if (processAllTypes || contains(simplexTypesToBenchmark, "dual"))
    runDualSimplexWithImplicitBounds<
        T, SimplexTraits<T, MatrixRepresentationType::NORMAL>>(
        *linearProgram, simplexTableauType, lpOptStatisticsVec);

  if (processAllTypes || contains(simplexTypesToBenchmark, "dual_sparse"))
    runDualSimplexWithImplicitBounds<
        T, SimplexTraits<T, MatrixRepresentationType::SPARSE>>(
        *linearProgram, simplexTableauType, lpOptStatisticsVec);

  if (absl::GetFlag(FLAGS_run_gomory)) {
    runDualSimplexGomoryWithPrimalCuts<
        T, SimplexTraits<T, MatrixRepresentationType::SPARSE>>(
        *linearProgram, simplexTableauType, lpOptStatisticsVec);
    runDualSimplexGomoryWithPrimalCuts<
        T, SimplexTraits<T, MatrixRepresentationType::NORMAL>>(
        *linearProgram, simplexTableauType, lpOptStatisticsVec);
  }
}

template <typename T>
LPOptStatistics<T>
readLPModelAndOptimizeByGurobi(const std::filesystem::path &modelFileMpsPath) {
  GurobiOptimizer gurobiOptimizer(absl::GetFlag(FLAGS_gurobi_log_file),
                                  modelFileMpsPath);
  return gurobiOptimizer.optimize<T>(LPOptimizationType::LINEAR_RELAXATION);
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
void printLPOptStats(const LPOptStatisticsVec<T> &lpOptStatisticsVec) {
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

  for (const SimplexTableauType simplexTableauType :
       absl::GetFlag(FLAGS_simplex_tableau_types)) {
    if (const auto lpModelFile = absl::GetFlag(FLAGS_lp_model_file);
        lpModelFile.has_value()) {
      readLPModelAndProcess<FloatingPointT>(std::filesystem::path{*lpModelFile},
                                            simplexTableauType,
                                            lpOptStatisticsVec);
      lpOptStatisticsVec.push_back(
          readLPModelAndOptimizeByGurobi<FloatingPointT>(
              std::filesystem::path{*lpModelFile}));
    } else if (const auto lpModelsDirectory =
                   absl::GetFlag(FLAGS_lp_models_directory);
               lpModelsDirectory.has_value())
      for (const auto &lpModelFileEntry :
           std::filesystem::directory_iterator(*lpModelsDirectory)) {
        readLPModelAndProcess<FloatingPointT>(
            lpModelFileEntry.path(), simplexTableauType, lpOptStatisticsVec);
        lpOptStatisticsVec.push_back(
            readLPModelAndOptimizeByGurobi<FloatingPointT>(
                lpModelFileEntry.path()));
      }
  }

  printLPOptStats(lpOptStatisticsVec);
}
