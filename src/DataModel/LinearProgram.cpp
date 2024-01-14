#include "src/DataModel/LinearProgram.h"

#include "src/Util/LPPrinter.h"
#include "src/Util/SpdlogHeader.h"

namespace {
template <typename T>
int countBoundsSpecified(const std::vector<std::optional<T>> &bounds) {
  return std::count_if(
      bounds.begin(), bounds.end(),
      [](const std::optional<T> &bound) { return bound.has_value(); });
}

template <typename T>
int countFreeVariables(const std::vector<std::optional<T>> &lowerBounds,
                       const std::vector<std::optional<T>> &upperBounds) {
  int freeVarCount = 0;
  for (int varIdx = 0; varIdx < lowerBounds.size(); ++varIdx) {
    if (!lowerBounds[varIdx].has_value() && !upperBounds[varIdx].has_value()) {
      ++freeVarCount;
    }
  }
  return freeVarCount;
}

int countIntegerVariables(const std::vector<VariableInfo> &variableInfos) {
  return std::count_if(variableInfos.begin(), variableInfos.end(),
                       [](const VariableInfo &variableInfo) {
                         return variableInfo._type == VariableType::INTEGER;
                       });
}

int countFixedVariables(const std::vector<VariableInfo> &variableInfos) {
  return std::count_if(
      variableInfos.begin(), variableInfos.end(),
      [](const VariableInfo &variableInfo) { return variableInfo._isFixed; });
}

const auto isInteger = [](const auto val) {
  double integerPart;
  return std::modf(val, &integerPart) == 0.0;
};
} // namespace

template <typename T> void LinearProgram<T>::convertToStandardForm() {
  int nonEqualityRows = 0;
  for (const auto &rowInfo : _rowInfos)
    if (!rowInfo.isEqualityRow())
      ++nonEqualityRows;

  const int variableCountAtTheStart = _variableInfos.size();
  const int newVariableCount = variableCountAtTheStart + nonEqualityRows;
  int newVariableIdx = variableCountAtTheStart;

  const auto newSlackLabel = [&]() {
    const std::string firstPattern = "S" + std::to_string(newVariableIdx);
    return (_variableLabelSet.find(firstPattern) == _variableLabelSet.end())
               ? firstPattern
               : firstPattern + Constants::SLACK_SUFFIX;
  };

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    _constraintMatrix[rowIdx].resize(newVariableCount);
    if (_rowInfos[rowIdx].isEqualityRow())
      continue;

    switch (_rowInfos[rowIdx]._type) {
    // TODO - check variable type integer or continuous
    // TODO - fix slack labeling
    case RowType::LESS_THAN_OR_EQUAL: {
      _constraintMatrix[rowIdx][newVariableIdx] = 1;
      break;
    }
    case RowType::GREATER_THAN_OR_EQUAL: {
      _constraintMatrix[rowIdx][newVariableIdx] = -1;
      break;
    }
    default:
      break;
    }
    const auto newSlackLabelStr = newSlackLabel();
    _variableInfos.push_back(VariableInfo{
        ._label = newSlackLabelStr,
        ._type = areAllCoefficientsInteger(rowIdx) ? VariableType::INTEGER
                                                   : VariableType::CONTINUOUS,
        ._isSlack = true});
    _variableLabelSet.insert(newSlackLabelStr);
    _variableLowerBounds.push_back(0.0);
    _variableUpperBounds.push_back(std::nullopt);
    _objective.push_back(0.0);
    ++newVariableIdx;
    _rowInfos[rowIdx]._type = RowType::EQUALITY;
  }
}

template <typename T>
bool LinearProgram<T>::areAllCoefficientsInteger(const int rowIdx) const {
  return isInteger(_rightHandSides[rowIdx]) &&
         std::all_of(_constraintMatrix[rowIdx].begin(),
                     _constraintMatrix[rowIdx].end(), isInteger) &&
         std::all_of(_variableInfos.begin(), _variableInfos.end(),
                     [&](const auto &variableInfo) {
                       const int varIdx = &variableInfo - &_variableInfos[0];
                       return _constraintMatrix[rowIdx][varIdx] == 0.0 ||
                              variableInfo._type == VariableType::INTEGER;
                     });
}

template <typename T>
std::optional<LinearProgram<T>>
LinearProgram<T>::dualProgram(const bool addObjValueFirstVar) const {
  LinearProgram<T> dualProgram;
  dualProgram._name = _name + "_DUAL";
  const int indexOffset = (int)addObjValueFirstVar;
  const size_t dualVarCount =
      indexOffset + _rightHandSides.size() + 2 * _objective.size();
  const size_t dualRowCount = indexOffset + _objective.size();
  const auto getDualVarIdV = [&](const int primalVarIdx) {
    return indexOffset + _rightHandSides.size() + primalVarIdx;
  };
  const auto getDualVarIdW = [&](const int primalVarIdx) {
    return indexOffset + _rightHandSides.size() + _objective.size() +
           primalVarIdx;
  };
  if (addObjValueFirstVar) {
    dualProgram._rightHandSides.push_back(0.0);
  }
  dualProgram._rightHandSides.insert(dualProgram._rightHandSides.end(),
                                     _objective.begin(), _objective.end());

  if (addObjValueFirstVar) {
    dualProgram._objective.push_back(0.0);
  }
  dualProgram._objective.insert(dualProgram._objective.end(),
                                _rightHandSides.begin(), _rightHandSides.end());
  dualProgram._objective.resize(dualVarCount);
  for (int primalVarIdx = 0; primalVarIdx < _variableInfos.size();
       ++primalVarIdx) {
    if (_variableLowerBounds[primalVarIdx].has_value()) {
      dualProgram._objective[getDualVarIdV(primalVarIdx)] =
          *_variableLowerBounds[primalVarIdx];
    }
    if (_variableUpperBounds[primalVarIdx].has_value()) {
      dualProgram._objective[getDualVarIdW(primalVarIdx)] =
          -*_variableUpperBounds[primalVarIdx];
    }
  }

  for (int dualVarIdx = 0; dualVarIdx < dualVarCount; ++dualVarIdx) {
    dualProgram._objective[dualVarIdx] = -dualProgram._objective[dualVarIdx];
  }

  dualProgram._constraintMatrix = transpose(_constraintMatrix);
  if (addObjValueFirstVar) {
    dualProgram._constraintMatrix.insert(dualProgram._constraintMatrix.begin(),
                                         std::vector<T>());
    dualProgram._constraintMatrix[0].push_back(1.0);
    dualProgram._constraintMatrix[0].insert(
        dualProgram._constraintMatrix[0].end(),
        dualProgram._objective.begin() + 1, dualProgram._objective.end());
    for (int dualVarIdx = 1; dualVarIdx < dualVarCount; ++dualVarIdx) {
      dualProgram._constraintMatrix[0][dualVarIdx] =
          -dualProgram._constraintMatrix[0][dualVarIdx];
    }
  }

  for (int rowIdx = indexOffset; rowIdx < dualProgram._constraintMatrix.size();
       ++rowIdx) {
    auto &currentRow = dualProgram._constraintMatrix[rowIdx];
    if (addObjValueFirstVar) {
      currentRow.insert(currentRow.begin(), 0.0);
    }
    currentRow.resize(dualVarCount);
    currentRow[getDualVarIdV(rowIdx - indexOffset)] = 1.0;
    currentRow[getDualVarIdW(rowIdx - indexOffset)] = -1.0;
  }

  dualProgram._variableInfos.resize(dualVarCount);

  const auto dualVarLabel = [&](const int dualVarIdx) -> std::string {
    if (addObjValueFirstVar && dualVarIdx == 0) {
      return "DUAL_OBJ";
    }
    if (dualVarIdx < indexOffset + _rightHandSides.size()) {
      return fmt::format("PI_{}", dualVarIdx - indexOffset);
    }
    if (dualVarIdx < indexOffset + _rightHandSides.size() + _objective.size()) {
      return fmt::format("V_{}",
                         dualVarIdx - _rightHandSides.size() - indexOffset);
    }

    return fmt::format("W_{}", dualVarIdx - _rightHandSides.size() -
                                   _objective.size() - indexOffset);
  };

  for (int dualVarIdx = 0; dualVarIdx < dualVarCount; ++dualVarIdx) {
    const bool isDualVarFree =
        dualVarIdx < _rightHandSides.size() + indexOffset;
    dualProgram._variableInfos[dualVarIdx] = VariableInfo{
        ._label = dualVarLabel(dualVarIdx),
        ._type = VariableType::CONTINUOUS,
        ._isObjectiveVar = (addObjValueFirstVar && dualVarIdx == 0),
        ._isFree = isDualVarFree};
    dualProgram._variableInfos[dualVarIdx]._isFree = isDualVarFree;
  }

  dualProgram._variableLowerBounds.resize(dualVarCount);
  dualProgram._variableUpperBounds.resize(dualVarCount);
  for (int primalVarIdx = 0; primalVarIdx < _variableInfos.size();
       ++primalVarIdx) {
    const int dualVarIdxV = getDualVarIdV(primalVarIdx);
    const int dualVarIdxW = getDualVarIdW(primalVarIdx);
    dualProgram._variableLowerBounds[dualVarIdxV] =
        dualProgram._variableLowerBounds[dualVarIdxW] = 0.0;

    if (!_variableLowerBounds[primalVarIdx].has_value()) {
      dualProgram._variableUpperBounds[dualVarIdxV] = 0.0;
      dualProgram._variableInfos[dualVarIdxV]._isFixed = true;
    }

    if (!_variableUpperBounds[primalVarIdx].has_value()) {
      dualProgram._variableUpperBounds[dualVarIdxW] = 0.0;
      dualProgram._variableInfos[dualVarIdxW]._isFixed = true;
    }
  }

  dualProgram._rowInfos.resize(dualRowCount);
  for (int dualRowIdx = 0; dualRowIdx < dualProgram._rowInfos.size();
       ++dualRowIdx) {
    dualProgram._rowInfos[dualRowIdx]._type = RowType::EQUALITY;
  }

  dualProgram.logGeneralInformation();
  return dualProgram;
}

template <typename T>
size_t LinearProgram<T>::getOriginalVariablesCount() const {
  const auto firstSlackVarIt =
      std::find_if(_variableInfos.begin(), _variableInfos.end(),
                   [](const auto &variableInfo) {
                     return variableInfo._isSlack || variableInfo._isArtificial;
                   });

  return firstSlackVarIt - _variableInfos.begin();
}

template <typename T> std::string LinearProgram<T>::toString() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);

  lpPrinter.printLineBreak();
  lpPrinter.printVariableInfos(std::nullopt);
  lpPrinter.printLineBreak();
  lpPrinter.printMatrixWithRHS({}, _constraintMatrix, _rightHandSides);
  lpPrinter.printLineBreak();

  return basicInformationStr() + lpPrinter.toString();
}
template <typename T>
std::string LinearProgram<T>::toStringLpSolveFormat() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printInLpSolveFormat(_constraintMatrix, _objective, _rightHandSides,
                                 _variableLowerBounds, _variableUpperBounds);
  return lpPrinter.toString();
}

template <typename T>
std::string LinearProgram<T>::basicInformationStr() const {
  std::ostringstream oss;

  oss << "NAME: " << _name << '\n';
  oss << "VARIABLE COUNT: " << _variableInfos.size() << '\n';
  oss << "ROW COUNT: " << _rowInfos.size() << '\n';
  oss << "CONSTRAINT MATRIX DIMENSIONS: " << _constraintMatrix.size() << " x "
      << _constraintMatrix.front().size() << '\n';

  return oss.str();
}

template <typename T> void LinearProgram<T>::logGeneralInformation() const {
  SPDLOG_INFO("NAME: {}", _name);
  SPDLOG_INFO(
      "VARIABLE COUNT: {} ROW COUNT: {}, CONSTRAINT MATRIX DIMENSIONS: {} x {}",
      _variableInfos.size(), _rowInfos.size(), _constraintMatrix.size(),
      _constraintMatrix.front().size());
  SPDLOG_INFO("LOWER BOUNDS {}, UPPER BOUNDS {}, FREE VARIABLES {}, FIXED "
              "VARIABLES {}, INTEGER "
              "VARIABLES {}",
              countBoundsSpecified(_variableLowerBounds),
              countBoundsSpecified(_variableUpperBounds),
              countFreeVariables(_variableLowerBounds, _variableUpperBounds),
              countFixedVariables(_variableInfos),
              countIntegerVariables(_variableInfos));
}

template <typename T>
bool LinearProgram<T>::checkIfAllBoundsAreSpeficied() const {
  return (_variableLowerBounds.size() == _variableInfos.size()) &&
         std::all_of(_variableLowerBounds.begin() + 1,
                     _variableLowerBounds.end(),
                     [&](const std::optional<T> &lb) {
                       const int varIdx = &lb - &_variableLowerBounds[0];
                       return _variableInfos[varIdx]._isSlack || lb.has_value();
                     }) &&
         (_variableUpperBounds.size() == _variableInfos.size()) &&
         std::all_of(_variableUpperBounds.begin() + 1,
                     _variableUpperBounds.end(),
                     [&](const std::optional<T> &ub) {
                       const int varIdx = &ub - &_variableUpperBounds[0];
                       return _variableInfos[varIdx]._isSlack || ub.has_value();
                     });
}

template <typename T>
bool LinearProgram<T>::isPureIP(const bool checkObjVar) const {
  auto startIter = _variableInfos.begin() + (checkObjVar ? 0 : 1);
  return std::all_of(startIter, _variableInfos.end(),
                     [](const VariableInfo &varInfo) {
                       return varInfo._type == VariableType::INTEGER;
                     });
}

template <typename T>
bool LinearProgram<T>::allCoefficientsAreIntegers() const {
  return std::all_of(_objective.begin(), _objective.end(), isInteger) &&
         std::all_of(_rightHandSides.begin(), _rightHandSides.end(),
                     isInteger) &&
         std::all_of(_constraintMatrix.begin(), _constraintMatrix.end(),
                     [&](const auto &coeffRow) {
                       return std::all_of(coeffRow.begin(), coeffRow.end(),
                                          isInteger);
                     });
}

template <typename T>
bool LinearProgram<T>::allVariablesAreNonnegative() const {
  return std::all_of(_variableLowerBounds.begin() + 1,
                     _variableLowerBounds.end(), [](const auto lowerBound) {
                       return lowerBound.has_value() && *lowerBound >= 0.0;
                     });
}

template class LinearProgram<double>;
template class LinearProgram<long double>;