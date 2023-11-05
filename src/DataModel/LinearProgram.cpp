#include "src/DataModel/LinearProgram.h"

#include "src/Util/LPPrinter.h"

template <typename T> void LinearProgram<T>::convertToStandardForm() {
  int nonEqualityRows = 0;
  for (const auto &rowInfo : _rowInfos)
    if (!rowInfo.isEqualityRow())
      ++nonEqualityRows;

  const int variableCountAtTheStart = _variableInfos.size();
  const int newVariableCount = variableCountAtTheStart + nonEqualityRows;
  int newVariableIdx = variableCountAtTheStart;
  _isVariableFreeBitset.resize(newVariableCount);

  const auto newSlackLabel = [&]() {
    const std::string firstPattern = "S" + std::to_string(newVariableIdx + 1);
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
    _variableInfos.push_back(VariableInfo{._label = newSlackLabelStr,
                                          ._type = VariableType::CONTINUOUS,
                                          ._isSlack = true});
    _variableLabelSet.insert(newSlackLabelStr);
    _variableLowerBounds.push_back(0.0);
    _variableUpperBounds.push_back(std::nullopt);
    ++newVariableIdx;
    _rowInfos[rowIdx]._type = RowType::EQUALITY;
  }
}

template <typename T>
size_t LinearProgram<T>::getOriginalVariablesCount() const {
  const auto firstSlackVarIt = std::find_if(
      _variableInfos.begin(), _variableInfos.end(),
      [](const auto &variableInfo) { return variableInfo._isSlack; });

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
  const auto isInteger = [](const auto val) {
    double integerPart;
    return std::modf(val, &integerPart) == 0.0;
  };

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