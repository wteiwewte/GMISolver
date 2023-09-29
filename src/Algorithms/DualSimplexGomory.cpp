#include "src/Algorithms/DualSimplexGomory.h"

#include "src/Algorithms/LexicographicOptimizer.h"
#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"

#include <fmt/format.h>

template <typename T, typename SimplexTraitsT>
DualSimplexGomory<T, SimplexTraitsT>::DualSimplexGomory(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const DualSimplexRowPivotRule dualSimplexRowPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _dualSimplexRowPivotRule(dualSimplexRowPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption) {}

template <typename T, typename SimplexTraitsT>
std::string DualSimplexGomory<T, SimplexTraitsT>::type() const {
  return "DUAL SIMPLEX GOMORY WITH PRIMAL CUTS (" +
         std::string(SimplexTraitsT::useSparseRepresentationValue ? "SPARSE"
                                                                  : "NORMAL") +
         ')';
}

template <typename T, typename SimplexTraitsT>
IPOptStatistics<T> DualSimplexGomory<T, SimplexTraitsT>::run(
    const LexicographicReoptType lexicographicReoptType,
    const LPOptimizationType lpOptimizationType,
    const GomoryCutChoosingRule gomoryCutChoosingRule) {
  SPDLOG_INFO("GOMORY WITH {} LEXICOGRAPHIC REOPTIMIZATION",
              lexicographicReoptTypeToStr(lexicographicReoptType));
  IPOptStatistics<T> ipOptStatistics;
  int relaxationNo = 1;
  ipOptStatistics._lpRelaxationStats.emplace_back() =
      runImpl(relaxationNo, lexicographicReoptType);

  if (lpOptimizationType == LPOptimizationType::LINEAR_RELAXATION)
    return ipOptStatistics;

  while (true) {
    ++relaxationNo;
    if (relaxationNo > 5)
      break;
    SPDLOG_INFO("{}TH GOMORY ROUND", relaxationNo);
    checkIfNonBasicVarsAreIntegral();
    const auto fractionalBasisRows =
        collectFractionalBasisRowIndices(gomoryCutChoosingRule);
    SPDLOG_INFO("FOUND {} FRACTIONAL VARIABLES - ROW IDXS [{}]",
                fractionalBasisRows.size(),
                fmt::join(fractionalBasisRows, ", "));
    if (fractionalBasisRows.empty())
      break;

    addCutRows(relaxationNo, fractionalBasisRows);
    addSlackVars(relaxationNo, fractionalBasisRows);
    _simplexTableau.initMatrixRepresentations();
    if (!_reinversionManager.reinverse())
      break;

    SPDLOG_INFO("AFTER REINVERSION");
    SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());
    SPDLOG_INFO(_simplexTableau.toStringSolution());

    ipOptStatistics._lpRelaxationStats.emplace_back() =
        runImpl(relaxationNo, lexicographicReoptType);
  }

  ipOptStatistics._optimalValue =
      ipOptStatistics._lpRelaxationStats.back()
          ._lexicographicReoptStats._objectiveValueAfterLexReopt;
  ipOptStatistics._optimalSolution = _simplexTableau._x;

  return ipOptStatistics;
}

template <typename T, typename SimplexTraitsT>
LPRelaxationStatistics<T> DualSimplexGomory<T, SimplexTraitsT>::runImpl(
    const int relaxationNo,
    const LexicographicReoptType lexicographicReoptType) {
  const auto relaxationId = [&relaxationNo] {
    return fmt::format("{}TH_RELAX", relaxationNo);
  };

  LPRelaxationStatistics<T> relaxationStats;
  relaxationStats._relaxationOptStats = dualSimplex().run(relaxationId());
  SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());
  //  SPDLOG_INFO(_simplexTableau.toStringSolution());

  relaxationStats._lexicographicReoptStats =
      lexicographicOptimizer().run(lexicographicReoptType, relaxationId());
  SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());
  //  SPDLOG_INFO(_simplexTableau.toStringSolution());
  return relaxationStats;
}

template <typename T, typename SimplexTraitsT>
RevisedDualSimplexPFIBounds<T, SimplexTraitsT>
DualSimplexGomory<T, SimplexTraitsT>::dualSimplex() const {
  return RevisedDualSimplexPFIBounds<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _dualSimplexRowPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
LexicographicOptimizer<T, SimplexTraitsT>
DualSimplexGomory<T, SimplexTraitsT>::lexicographicOptimizer() const {
  return LexicographicOptimizer<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _primalSimplexColumnPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomory<T, SimplexTraitsT>::checkIfNonBasicVarsAreIntegral()
    const {
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    const auto &varInfo = _simplexTableau._variableInfos[varIdx];
    if (varInfo._isArtificial || varInfo._isSlack || varInfo._isFixed)
      continue;

    if (varInfo._type == VariableType::INTEGER &&
        !NumericalTraitsT::isInteger(_simplexTableau._x[varIdx])) {
      if (!_simplexTableau._simplexBasisData
               ._isBasicColumnIndexBitset[varIdx]) {
        SPDLOG_WARN("NON-BASIC VAR IDX {} IS REQUIRED TO BE INTEGER BUT ITS "
                    "VALUE IS NOT INTEGER",
                    varIdx);
      }
    }
  }
}

template <typename T, typename SimplexTraitsT>
std::vector<int>
DualSimplexGomory<T, SimplexTraitsT>::collectFractionalBasisRowIndices(
    const GomoryCutChoosingRule gomoryCutChoosingRule) const {
  std::vector<int> fractionalBasisVarsRowIndices;
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    const auto basicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
    const auto &varInfo = _simplexTableau._variableInfos[basicVarIdx];
    // TODO maybe checking slack variables is needed
    // TODO for IP it doesn't seem so
    // TODO for MIP it does - to be verified
    if (varInfo._isArtificial || varInfo._isSlack || varInfo._isFixed)
      continue;

    if (varInfo._type == VariableType::INTEGER &&
        !NumericalTraitsT::isInteger(_simplexTableau._x[basicVarIdx])) {
      fractionalBasisVarsRowIndices.push_back(rowIdx);

      if (gomoryCutChoosingRule == GomoryCutChoosingRule::FIRST)
        return fractionalBasisVarsRowIndices;
    }
  }
  return fractionalBasisVarsRowIndices;
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomory<T, SimplexTraitsT>::addCutRows(
    const int relaxationNo,
    const std::vector<int> &fractionalBasisVarsRowIndices) const {
  const auto oldBasisSize = _simplexTableau._rowInfos.size();
  const auto newBasisSize = oldBasisSize + fractionalBasisVarsRowIndices.size();
  _simplexTableau._rowInfos.reserve(newBasisSize);
  _simplexTableau._constraintMatrix.reserve(newBasisSize);
  _simplexTableau._rightHandSides.reserve(newBasisSize);
  _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap.resize(
      newBasisSize);

  int rowNo = 0;
  for (const auto rowIdx : fractionalBasisVarsRowIndices) {
    _simplexTableau._rowInfos.push_back(
        RowInfo{._label = fmt::format("CUT_{}_{}", relaxationNo, rowNo),
                ._type = RowType::EQUALITY});

    const auto tableauRow = _simplexTableau.computeTableauRowGeneric(
        rowIdx); // FIXME this function probably doesnt work after adding new
                 // cuts
    auto &newCutRow = _simplexTableau._constraintMatrix.emplace_back();
    newCutRow.resize(_simplexTableau._variableInfos.size() +
                         fractionalBasisVarsRowIndices.size(),
                     0.0);
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (!_simplexTableau._simplexBasisData
               ._isBasicColumnIndexBitset[varIdx] &&
          !_simplexTableau._variableInfos[varIdx]._isArtificial) {
        newCutRow[varIdx] = std::floor(tableauRow[varIdx]) - tableauRow[varIdx];
      } else {
        //        SPDLOG_INFO("BASIC VAR IDX {} LABEL {}", varIdx,
        //                    _simplexTableau._variableInfos[varIdx]._label);
      }
      //      SPDLOG_INFO("CUT ROW - VAR IDX {} VAR LABEL {} VALUE {}", varIdx,
      //                  _simplexTableau._variableInfos[varIdx]._label,
      //                  newCutRow[varIdx]);
    }
    newCutRow[_simplexTableau._variableInfos.size() + rowNo] = 1.0;

    [[maybe_unused]] const auto basicColumnIdx =
        _simplexTableau.basicColumnIdx(rowIdx);
    //    SPDLOG_INFO("ROW IDX {} RHS VALUE {} SOL VALUE {}", rowIdx,
    //                _simplexTableau._rightHandSides[rowIdx],
    //                _simplexTableau._x[basicColumnIdx]);
    //    SPDLOG_INFO("RHS DIFF {}",
    //                std::floor(_simplexTableau._rightHandSides[rowIdx]) -
    //                    _simplexTableau._rightHandSides[rowIdx]);
    //    SPDLOG_INFO("SOL DIFF {}",
    //    std::floor(_simplexTableau._x[basicColumnIdx]) -
    //                                   _simplexTableau._x[basicColumnIdx]);
    _simplexTableau._rightHandSides.push_back(
        std::floor(_simplexTableau._rightHandSides[rowIdx]) -
        _simplexTableau._rightHandSides[rowIdx]);
    _simplexTableau._initialRightHandSides.push_back(
        std::floor(_simplexTableau._rightHandSides[rowIdx]) -
        _simplexTableau._rightHandSides[rowIdx]);
    _simplexTableau._simplexBasisData
        ._rowToBasisColumnIdxMap[oldBasisSize + rowNo] =
        _simplexTableau._variableInfos.size() + rowNo;
    ++rowNo;
  }
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomory<T, SimplexTraitsT>::addSlackVars(
    const int relaxationNo,
    const std::vector<int> &fractionalBasisVarsRowIndices) const {
  const auto oldVarCount = _simplexTableau._variableInfos.size();
  const auto newVarCount = oldVarCount + fractionalBasisVarsRowIndices.size();
  _simplexTableau._objectiveRow.resize(newVarCount, 0.0);
  _simplexTableau._reducedCosts.resize(newVarCount, 0.0);
  _simplexTableau._variableInfos.reserve(newVarCount);
  _simplexTableau._variableLowerBounds.reserve(newVarCount);
  _simplexTableau._variableUpperBounds.reserve(newVarCount);
  _simplexTableau._simplexBasisData.resizeVarCount(newVarCount);

  int newVarIdx = 0;
  const auto newSlackLabel = [&]() {
    const std::string firstPattern =
        fmt::format("CUT_{}_S_{}", relaxationNo, newVarIdx + 1);
    return (_simplexTableau._variableLabelSet.find(firstPattern) ==
            _simplexTableau._variableLabelSet.end())
               ? firstPattern
               : firstPattern + Constants::SLACK_SUFFIX;
  };

  for (; newVarIdx < fractionalBasisVarsRowIndices.size(); ++newVarIdx) {
    const auto newSlackLabelStr = newSlackLabel();
    _simplexTableau._variableInfos.push_back(VariableInfo{
        newSlackLabelStr, VariableType::INTEGER, true, false, false});
    _simplexTableau._variableLabelSet.insert(newSlackLabelStr);
    _simplexTableau._variableLowerBounds.push_back(0.0);
    _simplexTableau._variableUpperBounds.push_back(std::nullopt);
    _simplexTableau._simplexBasisData
        ._isBasicColumnIndexBitset[oldVarCount + newVarIdx] = true;
  }

  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    _simplexTableau._constraintMatrix[rowIdx].resize(newVarCount, 0.0);
  }
}

template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;