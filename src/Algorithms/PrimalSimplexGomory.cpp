#include "src/Algorithms/PrimalSimplexGomory.h"

#include "src/Algorithms/LexicographicOptimizer.h"
#include "src/Algorithms/PrimalSimplex.h"
#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexTableauResizer.h"
#include "src/Util/Time.h"

#include <random>

namespace {
std::vector<int> getRandomElement(const std::vector<int> &vec) {
  std::vector<int> out;
  std::sample(vec.begin(), vec.end(), std::back_inserter(out), 1,
              std::mt19937{std::random_device{}()});
  return out;
}
} // namespace

template <typename T, typename SimplexTraitsT>
PrimalSimplexGomory<T, SimplexTraitsT>::PrimalSimplexGomory(
    const LinearProgram<T> &primalLinearProgram,
    SimplexTableau<T, SimplexTraitsT> &dualSimplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption,
    const SlackCutRemovalCondition slackCutRemovalCondition,
    const LexicographicReoptType lexicographicReoptType,
    const int cutRoundLimit)
    : _primalLinearProgram(primalLinearProgram),
      _dualSimplexTableau(dualSimplexTableau),
      _reinversionManager(reinversionManager),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption),
      _slackCutRemovalCondition(slackCutRemovalCondition),
      _lexicographicReoptType(lexicographicReoptType),
      _cutRoundLimit(cutRoundLimit) {}

template <typename T, typename SimplexTraitsT>
std::string PrimalSimplexGomory<T, SimplexTraitsT>::type() const {
  return fmt::format(
      "PRIMAL SIMPLEX GOMORY WITH DUAL CUTS ({}, {})",
      simplexTableauTypeToStr(_dualSimplexTableau._simplexTableauType),
      matrixRepresentationTypeToStr(SimplexTraitsT::matrixRepresentationType));
}

template <typename T, typename SimplexTraitsT>
IPOptStatistics<T> PrimalSimplexGomory<T, SimplexTraitsT>::run(
    const LPOptimizationType lpOptimizationType,
    const GomoryCutChoosingRule gomoryCutChoosingRule) {
  SPDLOG_INFO("PRIMAL GOMORY WITH {} LEXICOGRAPHIC REOPTIMIZATION",
              lexicographicReoptTypeToStr(_lexicographicReoptType));
  IPOptStatistics<T> ipOptStatistics{
      ._lpName = _dualSimplexTableau.getName(),
      ._algorithmType = type(),
      ._reinversionFrequency = _reinversionManager.reinversionFrequency()};
  _dualSimplexTableau.setObjective(
      _dualSimplexTableau._initialProgram.getObjective());
  _dualSimplexTableau.calculateDualExplicit();
  //  if (_dualSimplexTableau._simplexTableauType == SimplexTableauType::FULL) {
  //    _dualSimplexTableau._basisMatrixInverse.clear();
  //  } TODO NEEDED FOR CUT GENERATION?

  ipOptStatistics._elapsedTimeSec = executeAndMeasureTime([&] {
    int relaxationNo = 1;
    ipOptStatistics._lpRelaxationStats.emplace_back() = runImpl(relaxationNo);

    if (lpOptimizationType == LPOptimizationType::LINEAR_RELAXATION)
      return;

    if (_dualSimplexTableau._simplexTableauType != SimplexTableauType::FULL) {
      SPDLOG_ERROR(
          "GOMORY CUTS NOT SUPPORTED FOR SIMPLEX TABLEAU TYPE {}",
          simplexTableauTypeToStr(_dualSimplexTableau._simplexTableauType));
      return;
    }

    while (true) {
      ++relaxationNo;
      if (relaxationNo > _cutRoundLimit)
        break;
      SPDLOG_INFO(
          _dualSimplexTableau.toStringSolutionWithDual(_primalLinearProgram));
      SPDLOG_INFO(_dualSimplexTableau.toStringObjectiveValue());

      SPDLOG_INFO("{}TH PRIMAL GOMORY ROUND", relaxationNo);
      const auto fractionalDualCoordinates =
          collectFractionalDualCoordinates(gomoryCutChoosingRule);
      SPDLOG_INFO("FOUND {} FRACTIONAL VARIABLES - DUAL VAR IDXS [{}]",
                  fractionalDualCoordinates.size(),
                  fmt::join(fractionalDualCoordinates, ", "));
      if (fractionalDualCoordinates.empty())
        break;

      addCutColumns(relaxationNo, fractionalDualCoordinates);

      _dualSimplexTableau.calculateSolution();
      _dualSimplexTableau.calculateCurrentObjectiveValue();

      SPDLOG_INFO("AFTER ADDITION OF NEW CUTS");
      SPDLOG_INFO(
          _dualSimplexTableau.toStringSolutionWithDual(_primalLinearProgram));
      SPDLOG_INFO(_dualSimplexTableau.toString());
      SPDLOG_INFO(_dualSimplexTableau.toStringObjectiveValue());

      ipOptStatistics._lpRelaxationStats.emplace_back() = runImpl(relaxationNo);
      SPDLOG_INFO("AFTER REOPT");
      SPDLOG_INFO(
          _dualSimplexTableau.toStringSolutionWithDual(_primalLinearProgram));
      SPDLOG_INFO(_dualSimplexTableau.toString());
      SPDLOG_INFO(_dualSimplexTableau.toStringObjectiveValue());
    }
  });

  ipOptStatistics._optimalValue = ipOptStatistics._lpRelaxationStats.back()
                                      ._lexicographicReoptStats._optimalValue;
  ipOptStatistics._optimalSolution = _dualSimplexTableau._x;
  ipOptStatistics._optResult = _dualSimplexTableau._result;

  return ipOptStatistics;
}

template <typename T, typename SimplexTraitsT>
LPRelaxationStatistics<T>
PrimalSimplexGomory<T, SimplexTraitsT>::runImpl(const int relaxationNo) {
  const auto relaxationId = [&relaxationNo] {
    return fmt::format("{}TH_RELAX", relaxationNo);
  };

  LPRelaxationStatistics<T> relaxationStats;
  relaxationStats._relaxationOptStats = primalSimplex().run(
      relaxationId(), PrintSimplexOptSummary::YES, PrimalPhase::TWO);
  _dualSimplexTableau.calculateSolution();
  _dualSimplexTableau.calculateCurrentObjectiveValue();
  SPDLOG_DEBUG(_dualSimplexTableau.toStringObjectiveValue());
  SPDLOG_DEBUG(_dualSimplexTableau.toStringSolution());

  relaxationStats._lexicographicReoptStats =
      lexicographicOptimizer().run(relaxationId());
  SPDLOG_DEBUG(_dualSimplexTableau.toStringObjectiveValue());
  SPDLOG_DEBUG(_dualSimplexTableau.toStringSolution());
  return relaxationStats;
}

template <typename T, typename SimplexTraitsT>
PrimalSimplex<T, SimplexTraitsT>
PrimalSimplexGomory<T, SimplexTraitsT>::primalSimplex() const {
  return PrimalSimplex<T, SimplexTraitsT>(
      _dualSimplexTableau, _reinversionManager, _primalSimplexColumnPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
LexicographicOptimizer<T, SimplexTraitsT>
PrimalSimplexGomory<T, SimplexTraitsT>::lexicographicOptimizer() const {
  return LexicographicOptimizer<T, SimplexTraitsT>(
      _dualSimplexTableau, _reinversionManager, _primalSimplexColumnPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption,
      _lexicographicReoptType);
}

template <typename T, typename SimplexTraitsT>
std::vector<int>
PrimalSimplexGomory<T, SimplexTraitsT>::collectFractionalDualCoordinates(
    const GomoryCutChoosingRule gomoryCutChoosingRule) const {
  std::vector<int> fractionalDualCoordinates;
  for (int dualIdx = 0;
       dualIdx < _primalLinearProgram.getOriginalVariablesCount(); ++dualIdx) {
    // TODO - check varInfo._type == VariableType::INTEGER?
    // TODO checking fixed?
    if (!NumericalTraitsT::isInteger(_dualSimplexTableau._y[dualIdx])) {
      fractionalDualCoordinates.push_back(dualIdx);
    }
  }

  switch (gomoryCutChoosingRule) {
  case GomoryCutChoosingRule::FIRST:
    return !fractionalDualCoordinates.empty()
               ? std::vector<int>{fractionalDualCoordinates[0]}
               : std::vector<int>{};
  case GomoryCutChoosingRule::ALL:
    return fractionalDualCoordinates;
  case GomoryCutChoosingRule::RANDOM:
    return getRandomElement(fractionalDualCoordinates);
  }
}

template <typename T, typename SimplexTraitsT>
void PrimalSimplexGomory<T, SimplexTraitsT>::addCutColumns(
    const int relaxationNo,
    const std::vector<int> &fractionalDualCoordinates) const {
  for (const auto dualCoord : fractionalDualCoordinates) {
    SPDLOG_INFO("FRACTION DUAL VAR IDX {} VALUE {}", dualCoord,
                _dualSimplexTableau._y[dualCoord]);
    std::vector<T> ithColumnOfBasisInverse(
        _dualSimplexTableau._rowInfos.size());
    for (int k = 0; k < _dualSimplexTableau._rowInfos.size(); ++k) {
      ithColumnOfBasisInverse[k] =
          _dualSimplexTableau._basisMatrixInverse[k][dualCoord];
    }

    std::vector<T> rVec(_dualSimplexTableau._rowInfos.size());
    for (int k = 0; k < _dualSimplexTableau._rowInfos.size(); ++k) {
      rVec[k] = -std::floor(ithColumnOfBasisInverse[k]);
    }

    std::vector<T> bWithTilde(_dualSimplexTableau._rowInfos.size());
    for (int k = 0; k < _dualSimplexTableau._rowInfos.size(); ++k) {
      T dotProduct{0};
      for (int j = 0; j < _dualSimplexTableau._rowInfos.size(); ++j) {
        const auto basicColumnIdx =
            _dualSimplexTableau._simplexBasisData._rowToBasisColumnIdxMap[j];
        dotProduct +=
            _dualSimplexTableau._constraintMatrix[k][basicColumnIdx] * rVec[j];
      }
      bWithTilde[k] = dotProduct + ((k == dualCoord) ? 1 : 0);
    }

    T yBWithTildeProduct{0};
    for (int k = 0; k < _dualSimplexTableau._rowInfos.size(); ++k) {
      SPDLOG_INFO("PART OF PRODUCT {}",
                  _dualSimplexTableau._y[k] * bWithTilde[k]);
      yBWithTildeProduct += _dualSimplexTableau._y[k] * bWithTilde[k];
    }
    SPDLOG_INFO("HALOOO {} {}", yBWithTildeProduct,
                std::floor(yBWithTildeProduct) - yBWithTildeProduct);
    _dualSimplexTableau._objectiveRow.push_back(std::floor(yBWithTildeProduct));
    _dualSimplexTableau._reducedCosts.push_back(std::floor(yBWithTildeProduct) -
                                                yBWithTildeProduct);

    const auto newCutVarLabel = [&]() {
      const std::string firstPattern =
          fmt::format("CUT_ROUND_{}_S_{}", relaxationNo, dualCoord + 1);
      return (_dualSimplexTableau._variableLabelSet.find(firstPattern) ==
              _dualSimplexTableau._variableLabelSet.end())
                 ? firstPattern
                 : firstPattern + Constants::SLACK_SUFFIX;
    };

    const auto newCutVarLabelStr = newCutVarLabel();
    _dualSimplexTableau._variableInfos.push_back(
        VariableInfo{._label = newCutVarLabelStr,
                     ._type = VariableType::CONTINUOUS,
                     ._isSlack = true, // TODO ???
                     ._isCutVar = true});

    _dualSimplexTableau._variableLabelSet.insert(newCutVarLabelStr);
    _dualSimplexTableau._variableLowerBounds.push_back(0);            // TODO ??
    _dualSimplexTableau._variableUpperBounds.push_back(std::nullopt); // TODO ??
    _dualSimplexTableau._simplexBasisData.resizeVarCount(
        _dualSimplexTableau._variableInfos.size());
    _dualSimplexTableau._simplexBasisData
        ._isColumnAtLowerBoundBitset[_dualSimplexTableau._variableInfos.size() -
                                     1] = true;
    for (int rowIdx = 0; rowIdx < _dualSimplexTableau._rowInfos.size();
         ++rowIdx) {
      _dualSimplexTableau._fullTableau[rowIdx].push_back(bWithTilde[rowIdx]);
    }
  }
}

template class PrimalSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class PrimalSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;