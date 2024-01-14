#include "src/Algorithms/LexicographicOptimizer.h"

#include "src/Algorithms/RevisedPrimalSimplexPFIBounds.h"
#include "src/Algorithms/SimplexTableau.h"
#include "src/Util/SpdlogHeader.h"

template <typename T, typename SimplexTraitsT>
LexicographicOptimizer<T, SimplexTraitsT>::LexicographicOptimizer(
    SimplexTableau<T, SimplexTraitsT> &simplexTableau,
    ReinversionManager<T, SimplexTraitsT> &reinversionManager,
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule,
    const int32_t objValueLoggingFrequency,
    const ValidateSimplexOption validateSimplexOption)
    : _simplexTableau(simplexTableau), _reinversionManager(reinversionManager),
      _primalSimplexColumnPivotRule(primalSimplexColumnPivotRule),
      _objValueLoggingFrequency(objValueLoggingFrequency),
      _validateSimplexOption(validateSimplexOption) {}

template <typename T, typename SimplexTraitsT>
RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>
LexicographicOptimizer<T, SimplexTraitsT>::primalSimplex() const {
  return RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>(
      _simplexTableau, _reinversionManager, _primalSimplexColumnPivotRule,
      _objValueLoggingFrequency, _validateSimplexOption);
}

template <typename T, typename SimplexTraitsT>
LexReoptStatistics<T> LexicographicOptimizer<T, SimplexTraitsT>::run(
    const LexicographicReoptType lexicographicReoptType,
    const std::string &lexOptId, const bool saveSolution) {
  LexReoptStatistics<T> lexReoptStats{._optResult = _simplexTableau._result};
  const auto initialIsFixedBitset = getIsFixedBitset();
  int curVarIdxToBeOptimized = 0;
  int varsFixedCount = initialIsFixedBitset.count();
  int optimizedVarCount = 0;
  fixNonBasicVariables(varsFixedCount);
  while (curVarIdxToBeOptimized < _simplexTableau._variableInfos.size() &&
         varsFixedCount < _simplexTableau._variableInfos.size()) {
    if (!_simplexTableau._variableInfos[curVarIdxToBeOptimized]._isFixed) {
      _simplexTableau.setObjective(
          singleVarObjective(curVarIdxToBeOptimized, lexicographicReoptType));
      auto lpStatisticsFromSingleVarOpt = primalSimplex().runImpl(
          fmt::format("{}_VAR_{}_{}", lexOptId, curVarIdxToBeOptimized,
                      lexicographicReoptTypeToStr(lexicographicReoptType)),
          false);
      lexReoptStats._lexLPReoptStatsVec.push_back(
          std::move(lpStatisticsFromSingleVarOpt));
      fixNonBasicVariables(varsFixedCount);
      ++optimizedVarCount;
    }
    ++curVarIdxToBeOptimized;
  }
  SPDLOG_INFO("AFTER LEXICOGRAPHIC {} OPTIMIZATION - FIXED VAR COUNT {} ALL "
              "VAR COUNT {}",
              lexicographicReoptTypeToStr(lexicographicReoptType),
              varsFixedCount, _simplexTableau._variableInfos.size());
  _simplexTableau.setObjective(_simplexTableau._initialProgram.getObjective());
  lexReoptStats._objectiveValueAfterLexReopt = _simplexTableau._objectiveValue;
  if (saveSolution) {
    const bool isFirstVarObjective =
        _simplexTableau._variableInfos[0]._isObjectiveVar;
    lexReoptStats._solution.assign(
        _simplexTableau._x.begin() + (isFirstVarObjective ? 1 : 0),
        _simplexTableau._x.begin() +
            _simplexTableau._initialProgram.getOriginalVariablesCount());
  }

  unfixVariables(initialIsFixedBitset);

  SPDLOG_INFO("LEXICOGRAPHIC {} REOPTIMIZATION RAN {} VARIABLE-SUBPROGRAMS "
              "(OUT OF {} ALL VARIABLES)",
              lexicographicReoptTypeToStr(lexicographicReoptType),
              optimizedVarCount, _simplexTableau._variableInfos.size());

  return lexReoptStats;
}

template <typename T, typename SimplexTraitsT>
boost::dynamic_bitset<>
LexicographicOptimizer<T, SimplexTraitsT>::getIsFixedBitset() const {
  boost::dynamic_bitset<> result(_simplexTableau._variableInfos.size());
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    result[varIdx] = _simplexTableau._variableInfos[varIdx]._isFixed;
  }
  return result;
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
void LexicographicOptimizer<T, SimplexTraitsT>::unfixVariables(
    const boost::dynamic_bitset<> &initialIsFixedBitset) {
  for (int varIdx = 0; varIdx < _simplexTableau._variableInfos.size();
       ++varIdx) {
    _simplexTableau._variableInfos[varIdx]._isFixed =
        initialIsFixedBitset[varIdx];
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
