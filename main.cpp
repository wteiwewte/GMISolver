#include "src/DataModel/SimplexTableu.h"
#include "src/Util/MpsReader.h"

//#include <absl/flags/flag.h>
//#include <absl/flags/parse.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <fmt/format.h>
#include <thread>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/ostream_sink.h"

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

//        spdlog::warn("HALOO");


    std::cout << "LP at the start\n";
    std::cout << linearProgram->toString() << '\n';
    std::cout << linearProgram->toStringLpSolveFormat() << '\n';

    linearProgram->convertToStandardForm();
    std::cout << "LP in standard form\n";
    std::cout << linearProgram->toString() << '\n';
    std::cout << linearProgram->toStringLpSolveFormat() << '\n';

    linearProgram->makeRightHandSidesNonNegative();
    std::cout << "LP in standard form with non-negative RHS\n";
    std::cout << linearProgram->toString() << '\n';
    std::cout << linearProgram->toStringLpSolveFormat() << '\n';

    linearProgram->addArtificialVariables();
    std::cout << "LP in standard form with artificial variables\n";
    std::cout << linearProgram->toString() << '\n';
    std::cout << linearProgram->toStringLpSolveFormat() << '\n';

    SimplexTableau<T> simplexTableau(*linearProgram, linearProgram->artificialObjectiveRow());
    std::cout << "Simplex tableau with artificial variables cost\n";
    std::cout << simplexTableau.toString() << '\n';
    std::cout << simplexTableau.toStringLpSolveFormat() << '\n';

    int iterCount = 1;
    while (true)
    {
        std::cout << "\nITERATION " << iterCount++ << "\n";
        const bool iterResult = simplexTableau.primalSimplex();
        std::cout << simplexTableau.toString() << '\n';

        if (iterResult)
            break;
    }
}

void setFileLogger()
{
    auto fileLogger = spdlog::basic_logger_mt("basic_logger", "logs.txt");
    spdlog::set_default_logger(fileLogger);
}

// TODO: comparing floating numbers with epsilon ?

int main(int argc, char** argv) {
//    const auto linearProgram = MpsReader::read<double>("/Users/janmelech/Desktop/test_1.mps");
//    const auto linearProgram = MpsReader::read<double>("/Users/janmelech/Desktop/miplib2017-testscript-v1.0.3/instances/miplib2017_ungzipped/seymour.mps");
//    absl::ParseCommandLine(argc, argv);
//    std::ostringstream oss;
//    auto ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_mt> (oss);
//    auto logger = std::make_shared<spdlog::logger>("my_logger", ostream_sink);
//    spdlog::set_default_logger(logger);
//    spdlog::info("HALOO");
    //    spdlog::shutdown();
    std::ofstream outFile("GMISolver_out.txt");
//    outFile << "";
    std::this_thread::sleep_for (std::chrono::seconds(1));
    std::cout.rdbuf(outFile.rdbuf());
    readLPModelAndProcess<double>("/Users/janmelech/Desktop/test_1.mps");
}
