#include "src/Algorithms/LexicographicOptimizer.h"

#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/SpdlogHeader.h"

template <typename T, typename SimplexTraitsT>
LexicographicOptimizer<T, SimplexTraitsT>::LexicographicOptimizer(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const int32_t objValueLoggingFrequency, const int32_t reinversionFrequency,
    const ValidateSimplex validateSimplex)
    : _simplexTableau(simplexTableau),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _reinversionFrequency(reinversionFrequency),
      _validateSimplex(validateSimplex) {}

template <typename T, typename SimplexTraitsT>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>
LexicographicOptimizer<T, SimplexTraitsT>::primalSimplex() const {
  return RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>(
      _simplexTableau, _primalSimplexColumnPivotRule, _objValueLoggingFrequency,
      _reinversionFrequency, _validateSimplex);
}

template <typename T, typename SimplexTraitsT>
LexReoptStatistics<T> LexicographicOptimizer<T, SimplexTraitsT>::run(
    const LexicographicReoptType lexicographicReoptType,
    const std::string &lexOptId) {
  LexReoptStatistics<T> lexReoptStats{._optResult = _simplexTableau._result};
  int curVarIdxToBeOptimized = 0;
  int varsFixedCount = 0;
  int optimizedVarCount = 0;
  while (curVarIdxToBeOptimized < _simplexTableau._variableInfos.size() &&
         varsFixedCount < _simplexTableau._variableInfos.size()) {
    fixNonBasicVariables(varsFixedCount);
    if (!_simplexTableau._variableInfos[curVarIdxToBeOptimized]._isFixed) {
      _simplexTableau.setObjective(
          singleVarObjective(curVarIdxToBeOptimized, lexicographicReoptType));
      auto lpStatisticsFromSingleVarOpt = primalSimplex().runImpl(
          fmt::format("{}_VAR_{}_{}", lexOptId, curVarIdxToBeOptimized,
                      (lexicographicReoptType == LexicographicReoptType::MIN
                           ? "MIN"
                           : "MAX")),
          false);
      lexReoptStats._lexLPReoptStatsVec.push_back(
          std::move(lpStatisticsFromSingleVarOpt));
      ++optimizedVarCount;
    }
    ++curVarIdxToBeOptimized;
  }
  _simplexTableau.setObjective(_simplexTableau._initialProgram.getObjective());
  lexReoptStats._objectiveValueAfterLexReopt = _simplexTableau._objectiveValue;
  unfixAllVariables();

  SPDLOG_INFO("LEXICOGRAPHIC {} REOPTIMIZATION RAN {} VARIABLE-SUBPROGRAMS",
              lexicographicReoptTypeToStr(lexicographicReoptType),
              optimizedVarCount);

  return lexReoptStats;
}
template <typename T, typename SimplexTraitsT>
std::vector<T> LexicographicOptimizer<T, SimplexTraitsT>::singleVarObjective(
    const int varIdx, const LexicographicReoptType lexicographicReoptType) {
  std::vector<T> result(_simplexTableau._variableInfos.size());
  result[varIdx] =
      lexicographicReoptType == LexicographicReoptType::MIN ? 1.0 : -1.0;
  return result;
}
template <typename T, typename SimplexTraitsT>
void LexicographicOptimizer<T, SimplexTraitsT>::fixNonBasicVariables(
    int &varsFixedCount) {
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (!_simplexTableau._variableInfos[varIdx]._isFixed &&
        !_simplexTableau._simplexBasisData._isBasicColumnIndexBitset[varIdx] &&
        !NumericalTraitsT::equal(_simplexTableau._reducedCosts[varIdx], 0.0)) {
      _simplexTableau._variableInfos[varIdx]._isFixed = true;
      ++varsFixedCount;
    }
  }
}
template <typename T, typename SimplexTraitsT>
void LexicographicOptimizer<T, SimplexTraitsT>::unfixAllVariables() {
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    if (_simplexTableau._variableLowerBounds[varIdx].has_value() &&
        _simplexTableau._variableLowerBounds[varIdx] ==
            _simplexTableau._variableUpperBounds[varIdx])
      continue;

    _simplexTableau._variableInfos[varIdx]._isFixed =
        false; // FIXME this unfixes fixed vars?
  }
}

template class LexicographicOptimizer<
    double, SimplexTraits<double, MatrixRepresentationType::SPARSE>>;
template class LexicographicOptimizer<
    double, SimplexTraits<double, MatrixRepresentationType::NORMAL>>;
template class LexicographicOptimizer<
    long double, SimplexTraits<long double, MatrixRepresentationType::SPARSE>>;
template class LexicographicOptimizer<
    long double, SimplexTraits<long double, MatrixRepresentationType::NORMAL>>;
