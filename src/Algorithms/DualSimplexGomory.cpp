#include "src/Algorithms/DualSimplexGomory.h"

#include "src/Algorithms/LexicographicOptimizer.h"
#include "src/Algorithms/ReinversionManager.h"
#include "src/Algorithms/RevisedDualSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Algorithms/SimplexTableauResizer.h"
#include "src/Util/Time.h"

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
  return fmt::format(
      "DUAL SIMPLEX GOMORY WITH PRIMAL CUTS ({}, {})",
      simplexTableauTypeToStr(_simplexTableau._simplexTableauType),
      matrixRepresentationTypeToStr(SimplexTraitsT::matrixRepresentationType));
}

template <typename T, typename SimplexTraitsT>
IPOptStatistics<T> DualSimplexGomory<T, SimplexTraitsT>::run(
    const LexicographicReoptType lexicographicReoptType,
    const LPOptimizationType lpOptimizationType,
    const GomoryCutChoosingRule gomoryCutChoosingRule) {
  SPDLOG_INFO("GOMORY WITH {} LEXICOGRAPHIC REOPTIMIZATION",
              lexicographicReoptTypeToStr(lexicographicReoptType));
  IPOptStatistics<T> ipOptStatistics{
      ._lpName = _simplexTableau.getName(),
      ._algorithmType = type(),
      ._reinversionFrequency = _reinversionManager.reinversionFrequency()};

  ipOptStatistics._elapsedTimeSec = executeAndMeasureTime([&] {
    int relaxationNo = 1;
    ipOptStatistics._lpRelaxationStats.emplace_back() =
        runImpl(relaxationNo, lexicographicReoptType);

    if (lpOptimizationType == LPOptimizationType::LINEAR_RELAXATION)
      return;

    if (_simplexTableau._simplexTableauType != SimplexTableauType::FULL) {
      SPDLOG_ERROR(
          "GOMORY CUTS NOT SUPPORTED FOR SIMPLEX TABLEAU TYPE {}",
          simplexTableauTypeToStr(_simplexTableau._simplexTableauType));
      return;
    }

    while (true) {
      ++relaxationNo;
      if (relaxationNo > 50)
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

      _simplexTableau.calculateCurrentObjectiveValue();
      _simplexTableau.calculateSolution();

      SPDLOG_INFO("AFTER ADDITION OF NEW CUTS");
      SPDLOG_INFO(_simplexTableau.toString());
      SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());
      SPDLOG_INFO(_simplexTableau.toStringSolution());

      ipOptStatistics._lpRelaxationStats.emplace_back() =
          runImpl(relaxationNo, lexicographicReoptType);

      if (removeCutsInBasis()) {
        SPDLOG_INFO("AFTER CUT REMOVAL");
        SPDLOG_INFO(_simplexTableau.toString());
        SPDLOG_INFO(_simplexTableau.toStringObjectiveValue());
        SPDLOG_INFO(_simplexTableau.toStringSolution());
        ipOptStatistics._lpRelaxationStats.emplace_back() =
            runImpl(relaxationNo, lexicographicReoptType);
      }
    }
  });

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
  SPDLOG_INFO(_simplexTableau.toStringSolution());
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
    const bool isOriginalVariable =
        basicVarIdx <
        _simplexTableau._initialProgram.getOriginalVariablesCount();
    if (isOriginalVariable && varInfo._type == VariableType::INTEGER &&
        !NumericalTraitsT::isInteger(_simplexTableau._x[basicVarIdx])) {
      if (varInfo._isFixed) {
        SPDLOG_ERROR("FIXED INTEGER VARIABLE {} HAS NON-INTEGER VALUE {}",
                     basicVarIdx, _simplexTableau._x[basicVarIdx]);
        continue;
      }

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
  _simplexTableau._fullTableau.reserve(newBasisSize);
  _simplexTableau._rightHandSides.reserve(newBasisSize);
  _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap.resize(
      newBasisSize);

  std::vector<VectorT> tableauRows;
  for (const auto rowIdx : fractionalBasisVarsRowIndices) {
    tableauRows.push_back(*_simplexTableau.computeTableauRowGeneric(rowIdx));
  }

  int rowNo = 0;
  for (const auto rowIdx : fractionalBasisVarsRowIndices) {
    _simplexTableau._rowInfos.push_back(
        RowInfo{._label = fmt::format("CUT_{}_{}", relaxationNo, rowNo),
                ._type = RowType::EQUALITY});

    const auto &tableauRow = tableauRows[rowNo];
    auto &newCutRow = _simplexTableau._fullTableau.emplace_back();
    newCutRow.resize(_simplexTableau._variableInfos.size() +
                     fractionalBasisVarsRowIndices.size());
    for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
         ++varIdx) {
      if (!_simplexTableau._simplexBasisData
               ._isBasicColumnIndexBitset[varIdx]) {
        if (_simplexTableau._simplexBasisData
                ._isColumnAtLowerBoundBitset[varIdx] ||
            _simplexTableau._isVariableFreeBitset[varIdx]) {
          newCutRow[varIdx] =
              std::floor(tableauRow[varIdx]) - tableauRow[varIdx];
        } else {
          SPDLOG_INFO("HALLOOO");
          newCutRow[varIdx] =
              -(tableauRow[varIdx] + std::floor(-tableauRow[varIdx]));
        }
      }
    }
    const int newCutVarIdx = _simplexTableau._variableInfos.size() + rowNo;
    newCutRow[newCutVarIdx] = 1.0;

    const auto cutRhs = std::floor(_simplexTableau._rightHandSides[rowIdx]) -
                        _simplexTableau._rightHandSides[rowIdx];
    _simplexTableau._rightHandSides.push_back(cutRhs);
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
  _simplexTableau._isVariableFreeBitset.resize(newVarCount);
  _simplexTableau._variableLowerBounds.reserve(newVarCount);
  _simplexTableau._variableUpperBounds.reserve(newVarCount);
  _simplexTableau._simplexBasisData.resizeVarCount(newVarCount);

  int newVarIdx = 0;
  const auto newSlackLabel = [&]() {
    const std::string firstPattern =
        fmt::format("CUT_ROUND_{}_S_{}", relaxationNo, newVarIdx + 1);
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
    _simplexTableau._fullTableau[rowIdx].resize(newVarCount, 0.0);
  }
}
template <typename T, typename SimplexTraitsT>
bool DualSimplexGomory<T, SimplexTraitsT>::removeCutsInBasis() const {
  bool atleastOneCutRemoved = false;
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    {
      const auto basicVarIdx = _simplexTableau.basicColumnIdx(rowIdx);
      const auto cutRowIdx =
          _simplexTableau._variableInfos[basicVarIdx]._cutRowIdx;
      // FIXME check if slack value == 0 ?
      //        if (NumericalTraitsT::greater(_simplexTableau._x[basicVarIdx],
      //        0.0))
      if (cutRowIdx.has_value()) {
        SPDLOG_INFO("FOUND BASIC CUT VAR IDX {}", basicVarIdx);
        removeCutFromBasis(rowIdx, basicVarIdx, *cutRowIdx);
        atleastOneCutRemoved = true;
      }
    }
  }
  return atleastOneCutRemoved;
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomory<T, SimplexTraitsT>::removeCutFromBasis(
    const int basisRowIdxMappedToCutVar, const int cutVarIdx,
    const int cutRowIdx) const {
  if (basisRowIdxMappedToCutVar != cutRowIdx) {
    //    SPDLOG_INFO(_simplexTableau.toString());
    //    SPDLOG_INFO(_simplexTableau.toStringSolution());
    dualSimplex().pivot(basisRowIdxMappedToCutVar, std::nullopt, true);
    dualSimplex().pivot(cutRowIdx, cutVarIdx, true);
    //    SPDLOG_INFO(_simplexTableau.toString());
    //    SPDLOG_INFO(_simplexTableau.toStringSolution());
    if (cutVarIdx != _simplexTableau.basicColumnIdx(cutRowIdx)) {
      SPDLOG_ERROR("CUT VAR IDX {} IS STILL NOT MAPPED TO CUR ROW IDX {}",
                   cutVarIdx, cutRowIdx);
    }
  }

  std::vector<bool> shouldVarBeRemoved(_simplexTableau._variableInfos.size(),
                                       false);
  std::vector<bool> shouldRowBeRemoved(_simplexTableau._rowInfos.size(), false);
  shouldVarBeRemoved[cutVarIdx] = true;
  shouldRowBeRemoved[cutRowIdx] = true;
  removeCuts(shouldVarBeRemoved, shouldRowBeRemoved);
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomory<T, SimplexTraitsT>::removeCuts(
    const std::vector<bool> &shouldVarBeRemoved,
    const std::vector<bool> &shouldRowBeRemoved) const {
  std::vector<int> oldRowIdxToNewRowIdx(_simplexTableau._rowInfos.size());
  int currentNewRowIdx = 0;
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    if (!shouldRowBeRemoved[rowIdx])
      oldRowIdxToNewRowIdx[rowIdx] = currentNewRowIdx++;
  }
  const auto newRowToBasicColumnIdxMap =
      _simplexTableau._simplexBasisData.fixMappingAfterRemoval(
          shouldVarBeRemoved, shouldRowBeRemoved);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_DEBUG
  debugLogOldAndNewBasis(oldRowIdxToNewRowIdx, newRowToBasicColumnIdxMap,
                         shouldRowBeRemoved);
#endif
  SimplexTableauResizer simplexTableauResizer(_simplexTableau,
                                              _reinversionManager);
  simplexTableauResizer.removeRows(shouldRowBeRemoved);
  simplexTableauResizer.removeVariables(shouldVarBeRemoved);

  if (_simplexTableau._simplexTableauType != SimplexTableauType::FULL) {
    _simplexTableau.initMatrixRepresentations();
  }
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (_simplexTableau._variableInfos[varIdx]._cutRowIdx.has_value())
      _simplexTableau._variableInfos[varIdx]._cutRowIdx =
          oldRowIdxToNewRowIdx[_simplexTableau._variableInfos[varIdx]
                                   ._cutRowIdx.value()];
  }
  _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap =
      newRowToBasicColumnIdxMap;
  _reinversionManager.reinverse();
}

template <typename T, typename SimplexTraitsT>
void DualSimplexGomory<T, SimplexTraitsT>::debugLogOldAndNewBasis(
    const std::vector<int> &oldRowIdxToNewRowIdx,
    const std::vector<int> &newRowToBasicColumnIdxMap,
    const std::vector<bool> &shouldRowBeRemoved) const {
  SPDLOG_DEBUG("OLD ROW->COLUMN MAPPING");
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    SPDLOG_DEBUG(
        "{}->{}", rowIdx,
        _simplexTableau._simplexBasisData._rowToBasisColumnIdxMap[rowIdx]);
  }
  SPDLOG_DEBUG("NEW ROW->COLUMN MAPPING");
  for (int rowIdx = 0; rowIdx < _simplexTableau._rowInfos.size(); ++rowIdx) {
    if (!shouldRowBeRemoved[rowIdx]) {
      [[maybe_unused]] const int newRowIdx = oldRowIdxToNewRowIdx[rowIdx];
      SPDLOG_DEBUG("{}->{}", newRowIdx, newRowToBasicColumnIdxMap[newRowIdx]);
    }
  }
}

template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class DualSimplexGomory<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;