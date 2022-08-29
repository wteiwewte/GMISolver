#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/PrimalSimplex.h"
#include "src/Util/MpsReader.h"

//#include <absl/flags/flag.h>
//#include <absl/flags/parse.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <spdlog/spdlog.h>
#include <thread>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/ostream_sink.h"
#include <fmt/format.h>

//ABSL_FLAG(std::string, model_file, "/Users/janmelech/Downloads/bm23.mps", "File with lp model`");

template <typename T>
void readLPModelAndProcess(const std::string& modelFileMps)
{
    auto linearProgram = MpsReader::read<T>(modelFileMps);
    if (!linearProgram.has_value())
    {
        spdlog::warn("Could not read properly lp from {} file", modelFileMps);
        return;
    }

    spdlog::info("LP at the start");
    spdlog::info(linearProgram->toString());
    spdlog::info(linearProgram->toStringLpSolveFormat());

    linearProgram->convertToStandardForm();
    spdlog::info("LP in standard form");
    spdlog::info(linearProgram->toString());
    spdlog::info(linearProgram->toStringLpSolveFormat());

    linearProgram->makeRightHandSidesNonNegative();
    spdlog::info("LP in standard form with non-negative RHS");
    spdlog::info(linearProgram->toString());
    spdlog::info(linearProgram->toStringLpSolveFormat());

    SimplexTableau<T> simplexTableau(*linearProgram);
    spdlog::info("Simplex tableau with artificial variables cost");
    spdlog::info(simplexTableau.toString());
    spdlog::info(simplexTableau.toStringLpSolveFormat());

    PrimalSimplex<T>(simplexTableau).runPhaseOne();
    spdlog::info("Simplex tableau after phase one");
    spdlog::info(simplexTableau.toString());
    spdlog::info(simplexTableau.toStringLpSolveFormat());
}

void initFileLogger()
{
    auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_st>("GMISolver_out.txt", true);
    auto fileLogger = std::make_shared<spdlog::logger>("fileLogger", fileSink);
    spdlog::set_default_logger(fileLogger);

    const std::string COUT_PATTERN = "%v";
    spdlog::set_pattern("%^[%Y-%m-%d %H:%M:%S.%e][%L]%$ %v");
    spdlog::set_level(spdlog::level::debug);
}

// TODO: comparing floating numbers with epsilon ?

int main(int argc, char** argv) {
    initFileLogger();
    std::this_thread::sleep_for (std::chrono::seconds(1));
    readLPModelAndProcess<double>("/Users/janmelech/Downloads/bm23.mps");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/test_1.mps");
}
