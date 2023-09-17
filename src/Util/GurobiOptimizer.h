#ifndef GMISOLVER_GUROBIOPTIMIZER_H
#define GMISOLVER_GUROBIOPTIMIZER_H

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
    if (_grbModel.get(GRB_IntAttr_IsMIP) != 0)
    {
      SPDLOG_INFO("The model {} is not a linear program", std::string{modelFileMpsPath});
    }
  }
  catch (GRBException e)
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

  void optimize()
  {
    try
    {
      _grbModel.optimize();

      int status = _grbModel.get(GRB_IntAttr_Status);

      if ((status == GRB_INF_OR_UNBD) || (status == GRB_INFEASIBLE) ||
          (status == GRB_UNBOUNDED))
      {
        SPDLOG_INFO("The model cannot be solved because it is infeasible or unbounded" );
        return;
      }

      if (status != GRB_OPTIMAL)
      {
        SPDLOG_INFO("Optimization was stopped with status {}", status );
        return;
      }

      double iterCount = _grbModel.get(GRB_DoubleAttr_IterCount);
      double runTime = _grbModel.get(GRB_DoubleAttr_Runtime);

      SPDLOG_INFO("Gurobi performed {} iterations in {} seconds", iterCount, runTime);
      SPDLOG_INFO("Model {} has optimal value {}", _grbModel.get(GRB_StringAttr_ModelName), _grbModel.get(GRB_DoubleAttr_ObjVal));
    }
    catch (GRBException e)
    {
      SPDLOG_INFO("Error code = {}", e.getErrorCode());
      SPDLOG_INFO(e.getMessage());
    }
    catch (...)
    {
      SPDLOG_INFO("Error during optimization");
    }
  }

private:
  GurobiEnvWrapper _grbEnvWrapper;
  GRBModel _grbModel;
};

#endif // GMISOLVER_GUROBIOPTIMIZER_H
