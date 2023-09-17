#ifndef GMISOLVER_GUROBIOPTIMIZER_H
#define GMISOLVER_GUROBIOPTIMIZER_H

#include "src/Util/LPOptStatistics.h"
#include "src/Util/SpdlogHeader.h"

#include <filesystem>

#include "gurobi_c++.h"

struct GurobiEnvWrapper
{
  GurobiEnvWrapper(const std::string& logFileName)
      : _grbEnv(true) {
    _grbEnv.set(GRB_StringParam_LogFile, logFileName);
    _grbEnv.set(GRB_IntParam_LogToConsole, 0);
    _grbEnv.start();
  }

  GRBEnv _grbEnv;
};

class GurobiOptimizer
{
public:
  GurobiOptimizer(const std::string& logFileName, const std::filesystem::path &modelFileMpsPath) try
      : _grbEnvWrapper(logFileName), _grbModel(_grbEnvWrapper._grbEnv, modelFileMpsPath)
  {
  }
  catch (const GRBException& e)
  {
    SPDLOG_INFO("Error code = {}", e.getErrorCode());
    SPDLOG_INFO(e.getMessage());
  }
  catch (...)
  {
    SPDLOG_INFO("Error during initialization");
  }

  void relaxModel()
  {
    auto vars = _grbModel.getVars();
    for (int varIdx = 0; varIdx < _grbModel.get(GRB_IntAttr_NumVars); ++varIdx)
    {
      vars[varIdx].set(GRB_CharAttr_VType, 'C');
    }
  }

  LPOptStatistics<double> optimize()
  {
    LPOptStatistics<double> lpOptStatistics{._lpName = _grbModel.get(GRB_StringAttr_ModelName), ._simplexAlgorithmType="GUROBI"};
    try
    {
      _grbModel.optimize();

      LPOptimizationResult gurobiLpOptResult = gurobiStatusToLPOptResult(_grbModel.get(GRB_IntAttr_Status));
      double iterCount = _grbModel.get(GRB_DoubleAttr_IterCount);
      double runTime = _grbModel.get(GRB_DoubleAttr_Runtime);
      double optimalValue = _grbModel.get(GRB_DoubleAttr_ObjVal);

      SPDLOG_INFO("Gurobi performed {} iterations in {} seconds", iterCount, runTime);
      SPDLOG_INFO("Model optimization status {}", lpOptimizationResultToStr(gurobiLpOptResult));
      SPDLOG_INFO("Model {} has optimal value {}", _grbModel.get(GRB_StringAttr_ModelName), optimalValue);

      lpOptStatistics._optResult = gurobiLpOptResult;
      lpOptStatistics._iterationCount = iterCount;
      lpOptStatistics._optimalValue = optimalValue;
    }
    catch (const GRBException& e)
    {
      SPDLOG_INFO("Error code = {}", e.getErrorCode());
      SPDLOG_INFO(e.getMessage());
    }
    catch (...)
    {
      SPDLOG_INFO("Error during optimization");
    }

    return lpOptStatistics;
  }

private:

  LPOptimizationResult gurobiStatusToLPOptResult(const int gurobiStatus)
  {
    switch(gurobiStatus) {
      case GRB_OPTIMAL:
        return LPOptimizationResult::BOUNDED_AND_FEASIBLE;
      case GRB_INFEASIBLE:
        return LPOptimizationResult::INFEASIBLE;
      case GRB_UNBOUNDED:
        return LPOptimizationResult::UNBOUNDED;
      case GRB_INF_OR_UNBD:
        return LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED;
      case GRB_ITERATION_LIMIT:
        return LPOptimizationResult::REACHED_ITERATION_LIMIT;
    }

    return LPOptimizationResult::UNKNOWN;
  }

  GurobiEnvWrapper _grbEnvWrapper;
  GRBModel _grbModel;
};

#endif // GMISOLVER_GUROBIOPTIMIZER_H
