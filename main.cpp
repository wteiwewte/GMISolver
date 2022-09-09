#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/PrimalSimplex.h"
#include "src/Algorithms/RevisedPrimalSimplexPFI.h"
#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "src/Util/MpsReader.h"

#include <iostream>
#include <filesystem>
#include <fstream>
#include <optional>
#include <vector>
#include <spdlog/spdlog.h>
#include <thread>

#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <spdlog/sinks/basic_file_sink.h>

ABSL_FLAG(std::string, log_level, "info", "log level");
ABSL_FLAG(std::string, log_file, "GMISolver_out.txt", "log file");
ABSL_FLAG(std::optional<std::string>, lp_model_file, std::nullopt, "Path to single lp model");
ABSL_FLAG(std::optional<std::string>, lp_models_directory, std::nullopt, "Path to directory with lp models");

template <typename T>
void runPrimalSimplexWithImplicitBounds(const LinearProgram<T>& linearProgram)
{
  SimplexTableau<T> simplexTableau(linearProgram, true);
  RevisedPrimalSimplexPFIBounds<T>(simplexTableau).run();
  SPDLOG_INFO(simplexTableau.toStringShortWithSolution());
}

template <typename T>
void runDualSimplexWithImplicitBounds(const LinearProgram<T>& linearProgram)
{
  SimplexTableau<T> simplexTableau(linearProgram, false);
  RevisedDualSimplexPFIBounds<T>(simplexTableau).run();
  SPDLOG_INFO(simplexTableau.toStringShortWithSolution());
}

template <typename T>
void readLPModelAndProcess(const std::string& modelFileMps)
{
  SPDLOG_INFO("Processing lp model from file {}", modelFileMps);
  auto linearProgram = MpsReader::read<T>(modelFileMps);
  if (!linearProgram.has_value())
  {
      spdlog::warn("Could not read properly lp from {} file", modelFileMps);
      return;
  }

  runPrimalSimplexWithImplicitBounds(*linearProgram);
  runDualSimplexWithImplicitBounds(*linearProgram);
}

void initFileLogger()
{
  auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_st>(absl::GetFlag(FLAGS_log_file), true);
  auto fileLogger = std::make_shared<spdlog::logger>("fileLogger", fileSink);
  spdlog::set_default_logger(fileLogger);

  const std::string COUT_PATTERN = "%v";
  spdlog::set_pattern("%^[%Y-%m-%d %H:%M:%S.%e][%L]%$ %v");
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_every(std::chrono::seconds(5));
}

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  initFileLogger();
  std::this_thread::sleep_for (std::chrono::seconds(1));

  if (const auto lpModelFile = absl::GetFlag(FLAGS_lp_model_file); lpModelFile.has_value())
    readLPModelAndProcess<double>(*lpModelFile);
  else if(const auto lpModelsDirectory = absl::GetFlag(FLAGS_lp_models_directory); lpModelsDirectory.has_value())
    for (const auto& lpModelFileEntry : std::filesystem::directory_iterator(*lpModelsDirectory))
      readLPModelAndProcess<double>(lpModelFileEntry.path());
}
