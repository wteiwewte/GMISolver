#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/PrimalSimplex.h"
#include "src/Algorithms/RevisedPrimalSimplexPFI.h"
#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"
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

//ABSL_FLAG(std::string, model_file, "/Users/janmelech/Desktop/magisterka/dane testowe/bm23.mps", "File with lp model`");

template <typename T>
void runPrimalSimplex(LinearProgram<T>& linearProgram)
{
  linearProgram.addBoundsToMatrix();
  linearProgram.convertToStandardForm();
  spdlog::info("LP in standard form");
  //    spdlog::info(linearProgram->toString());
  //    spdlog::info(linearProgram->toStringLpSolveFormat());

  linearProgram.makeRightHandSidesNonNegative();
  spdlog::info("LP in standard form with non-negative RHS");
  //    spdlog::info(linearProgram->toString());
  //    spdlog::info(linearProgram->toStringLpSolveFormat());

  SimplexTableau<T> simplexTableau(linearProgram, true);
  spdlog::info("Simplex tableau with artificial variables cost");
  //    spdlog::info(simplexTableau.toString());
  //    spdlog::info(simplexTableau.toStringLpSolveFormat());

  RevisedPrimalSimplexPFI<T>(simplexTableau).runPhaseOne();
  //    RevisedPrimalSimplexPFI<T>(simplexTableau).runPhaseOne();
  //    PrimalSimplex<T>(simplexTableau).runPhaseOne();
  spdlog::info("Simplex tableau after primal simplex");
  //    spdlog::info(simplexTableau.toString());
  spdlog::info(simplexTableau.toStringShortWithSolution());
}

template <typename T>
void runPrimalSimplexWithImplicitBounds(LinearProgram<T>& linearProgram)
{
  linearProgram.convertToStandardForm();
  spdlog::info("LP in standard form");
  //    spdlog::info(linearProgram->toString());
  //    spdlog::info(linearProgram->toStringLpSolveFormat());

  linearProgram.makeRightHandSidesNonNegative();
  spdlog::info("LP in standard form with non-negative RHS");
  //    spdlog::info(linearProgram->toString());
  //    spdlog::info(linearProgram->toStringLpSolveFormat());

  SimplexTableau<T> simplexTableau(linearProgram, true);
  spdlog::info("Simplex tableau with artificial variables cost");
  //    spdlog::info(simplexTableau.toString());
  //    spdlog::info(simplexTableau.toStringLpSolveFormat());

  RevisedPrimalSimplexPFIBounds<T>(simplexTableau).runPhaseOne();
  //    RevisedPrimalSimplexPFI<T>(simplexTableau).runPhaseOne();
  //    PrimalSimplex<T>(simplexTableau).runPhaseOne();
  spdlog::info("Simplex tableau after primal simplex");
  //    spdlog::info(simplexTableau.toString());
  spdlog::info(simplexTableau.toStringShortWithSolution());
}

template <typename T>
void runDualSimplexWithImplicitBounds(LinearProgram<T>& linearProgram)
{
  linearProgram.convertToStandardForm();
//  spdlog::info("LP in standard form");
//  //    spdlog::info(linearProgram->toString());
//  spdlog::info(linearProgram.toStringLpSolveFormat());
//
//  linearProgram.makeRightHandSidesNonNegative();
//  spdlog::info("LP in standard form with non-negative RHS");
//  //    spdlog::info(linearProgram->toString());
//  //    spdlog::info(linearProgram->toStringLpSolveFormat());
//
  SimplexTableau<T> simplexTableau(linearProgram, false);
//  spdlog::info("Simplex tableau with artificial variables cost");
//  //    spdlog::info(simplexTableau.toString());
//  spdlog::info(simplexTableau.toStringLpSolveFormat());
//
//  //    RevisedPrimalSimplexPFI<T>(simplexTableau).runPhaseOne();
//  //    PrimalSimplex<T>(simplexTableau).runPhaseOne();
//  spdlog::info("Simplex tableau after primal simplex");
//  //    spdlog::info(simplexTableau.toString());
  RevisedDualSimplexPFIBounds<T>(simplexTableau).run();
  spdlog::info(simplexTableau.toStringShortWithSolution());
}

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
    //    spdlog::info(linearProgram->toString());
    //    spdlog::info(linearProgram->toStringLpSolveFormat());

//    runPrimalSimplex(*linearProgram);
//    runPrimalSimplexWithImplicitBounds(*linearProgram);
    runDualSimplexWithImplicitBounds(*linearProgram);
}

void initFileLogger()
{
    auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_st>("GMISolver_out.txt", true);
    auto fileLogger = std::make_shared<spdlog::logger>("fileLogger", fileSink);
    spdlog::set_default_logger(fileLogger);

    const std::string COUT_PATTERN = "%v";
    spdlog::set_pattern("%^[%Y-%m-%d %H:%M:%S.%e][%L]%$ %v");
    spdlog::set_level(spdlog::level::info);
    spdlog::flush_every(std::chrono::seconds(5));
}

// TODO: comparing floating numbers with epsilon ?

int main(int argc, char** argv) {
    initFileLogger();
    std::this_thread::sleep_for (std::chrono::seconds(1));
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/test_1.mps");
    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/bm23");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/air01");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/air02");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/air03");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/air04");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/air05");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/diamond");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/bell5");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/fixnet4");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/sentoy");
//    readLPModelAndProcess<double>("/Users/janmelech/Desktop/magisterka/dane testowe/miplib2/miplib/mod010");
}
