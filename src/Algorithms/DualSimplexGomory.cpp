#include "src/Algorithms/DualSimplexGomory.h"

#include "src/Algorithms/LexicographicOptimizer.h"
#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexTableauResizer.h"

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
    if (relaxationNo > 25)
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
    SPDLOG_INFO(_simplexTableau.toString());
    SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());
    SPDLOG_INFO(_simplexTableau.toStringSolution());

    ipOptStatistics._lpRelaxationStats.emplace_back() =
        runImpl(relaxationNo, lexicographicReoptType);

    SPDLOG_INFO(_simplexTableau.toString());
    if (!removeSlackCuts())
      break;
    SPDLOG_INFO("AFTER CUT REMOVAL");
    SPDLOG_INFO(_simplexTableau.toString());
    SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());
    SPDLOG_INFO(_simplexTableau.toStringSolution());
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
  std::vector<VectorT> tableauRows;
  for (const auto rowIdx : fractionalBasisVarsRowIndices) {
    tableauRows.push_back(_simplexTableau.computeTableauRowGeneric(rowIdx));
  }

  for (const auto rowIdx : fractionalBasisVarsRowIndices) {
    _simplexTableau._rowInfos.push_back(
        RowInfo{._label = fmt::format("CUT_{}_{}", relaxationNo, rowNo),
                ._type = RowType::EQUALITY});

    //    const auto tableauRow = _simplexTableau.computeTableauRowGeneric(
    //        rowIdx); // FIXME this function probably doesnt work after adding
    //        new
    //                 // cuts
    const auto &tableauRow = tableauRows[rowNo];
    auto &newCutRow = _simplexTableau._constraintMatrix.emplace_back();
    newCutRow.resize(_simplexTableau._variableInfos.size() +
                         fractionalBasisVarsRowIndices.size(),
                     0.0);
    T rhsAddition = 0;
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (!_simplexTableau._simplexBasisData
               ._isBasicColumnIndexBitset[varIdx]) {
        const auto val = std::floor(tableauRow[varIdx]) - tableauRow[varIdx];
        if (_defineCutsInTermsOfOriginalVariables) {
          if (const auto cutRowIdx =
                  _simplexTableau._variableInfos[varIdx]._cutRowIdx;
              cutRowIdx.has_value()) {
            for (int varIdx2 = 0;
                 varIdx2 < _simplexTableau._variableInfos.size(); ++varIdx2) {
              if (varIdx2 != varIdx) {
                newCutRow[varIdx2] +=
                    (-val) *
                    _simplexTableau._constraintMatrix[*cutRowIdx][varIdx2];
              }
            }
            rhsAddition =
                (-val) * _simplexTableau._initialRightHandSides[*cutRowIdx];
          } else {
            newCutRow[varIdx] += val;
          }
        } else {
          newCutRow[varIdx] = val;
        }
      }
    }
    const int newCutVarIdx = _simplexTableau._variableInfos.size() + rowNo;
    newCutRow[newCutVarIdx] = 1.0;

    const auto cutRhs = rhsAddition +
                        std::floor(_simplexTableau._rightHandSides[rowIdx]) -
                        _simplexTableau._rightHandSides[rowIdx];
    _simplexTableau._rightHandSides.push_back(cutRhs);
    _simplexTableau._initialRightHandSides.push_back(cutRhs);
    _simplexTableau._simplexBasisData
        ._rowToBasisColumnIdxMap[oldBasisSize + rowNo] = newCutVarIdx;
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
    const int cutRowIdx = _simplexTableau._rowInfos.size() -
                          fractionalBasisVarsRowIndices.size() + newVarIdx;
    const auto newSlackLabelStr = newSlackLabel();
    _simplexTableau._variableInfos.push_back(
        VariableInfo{._label = newSlackLabelStr,
                     ._type = VariableType::INTEGER,
                     ._isSlack = true,
                     ._cutRowIdx = cutRowIdx});
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

template <typename T, typename SimplexTraitsT>
bool DualSimplexGomory<T, SimplexTraitsT>::removeSlackCuts() const {
  SPDLOG_INFO("INITIAL ROW COUNT {} CURRENT ROW COUNT {}",
              _simplexTableau._initialProgram.getRowInfos().size(),
              _simplexTableau._rowInfos.size());
  bool atleastOneRowShouldBeRemoved = false;
  std::vector<bool> shouldVarBeRemoved(_simplexTableau._variableInfos.size(),
                                       false);
  std::vector<bool> shouldRowBeRemoved(_simplexTableau._rowInfos.size(), false);
  for (int rowIdx = _simplexTableau._initialProgram.getRowInfos().size();
       rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    {
      const auto basicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
      const auto cutRowIdx =
          _simplexTableau._variableInfos[basicVarIdx]._cutRowIdx;
      if (cutRowIdx.has_value()) {
        // FIXME check if slack value == 0 ?
        //        if (NumericalTraitsT::greater(_simplexTableau._x[basicVarIdx],
        //        0.0))
        atleastOneRowShouldBeRemoved = true;
        shouldVarBeRemoved[basicVarIdx] = true;
        shouldRowBeRemoved[*cutRowIdx] = true;
      }
    }
  }
  if (!atleastOneRowShouldBeRemoved)
    return true;

  std::vector<int> oldRowIdxToNewRowIdx(_simplexTableau._rowInfos.size());
  int currentNewRowIdx = 0;
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    if (!shouldRowBeRemoved[rowIdx])
      oldRowIdxToNewRowIdx[rowIdx] = currentNewRowIdx++;
  }
  SPDLOG_INFO("OLD ROW->COLUMN MAPPING");
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    SPDLOG_INFO(
        "{}->{}", rowIdx,
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx]);
  }
  const auto newRowToBasicColumnIdxMap =
      _simplexTableau._simplexBasisData.fixMappingAfterRemoval(
          shouldVarBeRemoved, shouldRowBeRemoved);
  SPDLOG_INFO("NEW ROW->COLUMN MAPPING");
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    if (!shouldRowBeRemoved[rowIdx]) {
      const int newRowIdx = oldRowIdxToNewRowIdx[rowIdx];
      SPDLOG_INFO("{}->{}", newRowIdx, newRowToBasicColumnIdxMap[newRowIdx]);
    }
  }
  SimplexTableauResizer simplexTableauResizer(_simplexTableau,
                                              _reinversionManager);
  simplexTableauResizer.removeRows(shouldRowBeRemoved);
  simplexTableauResizer.removeVariables(shouldVarBeRemoved);
  _simplexTableau.initMatrixRepresentations();
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (_simplexTableau._variableInfos[varIdx]._cutRowIdx.has_value())
      _simplexTableau._variableInfos[varIdx]._cutRowIdx =
          oldRowIdxToNewRowIdx[_simplexTableau._variableInfos[varIdx]
                                   ._cutRowIdx.value()];
  }
  _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap =
      newRowToBasicColumnIdxMap;
  return _reinversionManager.reinverse();
}

template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;