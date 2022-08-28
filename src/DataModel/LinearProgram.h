#ifndef GMISOLVER_LINEARPROGRAM_H
#define GMISOLVER_LINEARPROGRAM_H

#include "CommonConstants.h"
#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/EnumTypes.h"
#include "src/Util/LPPrinter.h"

#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

template <typename T> class LinearProgram {
public:
  std::string toString() const;
  std::string toStringLpSolveFormat() const;

  void convertToStandardForm();
  void makeRightHandSidesNonNegative();
  void addArtificialVariables();

  size_t initialVariableCountInStandardForm() const {
    return _initialVariableCountInStandardForm;
  }

  void setInitialVariableCountInStandardForm(const int val) {
    _initialVariableCountInStandardForm = val;
  }

  std::vector<T> artificialObjectiveRow() const {
    std::vector<T> result(_constraintMatrix.front().size());

    for (int variableIdx = 0; variableIdx < _variableInfos.size();
         ++variableIdx)
      if (_variableInfos[variableIdx]._isArtificial)
        result[variableIdx] = 1;

    return result;
  }

private:
  friend struct MpsReader;

  template <typename U> friend class SimplexTableau;

  std::string basicInformationStr() const {
    std::ostringstream oss;

    oss << "NAME: " << _name << '\n';
    oss << "VARIABLE COUNT: " << _variableInfos.size() << '\n';
    oss << "ROW COUNT: " << _rowInfos.size() << '\n';
    oss << "CONSTRAINT MATRIX DIMENSIONS: " << _constraintMatrix.size() << " x "
        << _constraintMatrix.front().size() << '\n';

    return oss.str();
  }

  std::string _name;

  std::vector<RowInfo> _rowInfos;
  std::vector<VariableInfo> _variableInfos;
  std::set<std::string> _variableLabelSet;

  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;
  std::vector<T> _objectiveCoeffs;
  int _initialVariableCountInStandardForm{0};
};

template <typename T> std::string LinearProgram<T>::toString() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);

  lpPrinter.printLineBreak();
  lpPrinter.printVariableInfos();
  lpPrinter.printLineBreak();
  lpPrinter.printMatrixWithRHS({}, _constraintMatrix, _rightHandSides);
  lpPrinter.printLineBreak();

  return basicInformationStr() + lpPrinter.toString();
}
template <typename T>
std::string LinearProgram<T>::toStringLpSolveFormat() const {
  LPPrinter lpPrinter(_variableInfos, _rowInfos);
  lpPrinter.printInLpSolveFormat(_constraintMatrix, _rightHandSides);
  return lpPrinter.toString();
}

template <typename T> void LinearProgram<T>::convertToStandardForm() {

  int nonEqualityRows = 0;
  for (const auto &rowInfo : _rowInfos)
    if (!rowInfo.isEqualityRow() && !rowInfo.isObjectiveRow())
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
    if (_rowInfos[rowIdx].isEqualityRow() || _rowInfos[rowIdx].isObjectiveRow())
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
    _variableInfos.push_back(
        VariableInfo{newSlackLabelStr, VariableType::CONTINUOUS, true});
    _variableLabelSet.insert(newSlackLabelStr);
    ++newVariableIdx;
    _rowInfos[rowIdx]._type = RowType::EQUALITY;
  }

  setInitialVariableCountInStandardForm(_variableInfos.size());
}

template <typename T> void LinearProgram<T>::makeRightHandSidesNonNegative() {
  for (int rowIdx = 1; rowIdx < _rowInfos.size(); ++rowIdx) {
    if (_rightHandSides[rowIdx] < 0.0) {
      for (int variableIdx = 0; variableIdx < _variableInfos.size();
           ++variableIdx)
        _constraintMatrix[rowIdx][variableIdx] =
            -_constraintMatrix[rowIdx][variableIdx];
      _rightHandSides[rowIdx] = -_rightHandSides[rowIdx];
    }
  }
}

template <typename T> void LinearProgram<T>::addArtificialVariables() {
  const int variableCountAtTheStart = _variableInfos.size();
  const int newVariableCount = variableCountAtTheStart + _rowInfos.size() - 1;

  const auto newArtificialLabel = [&](const auto varIdx) {
    const std::string firstPattern = "A" + std::to_string(varIdx);
    return (_variableLabelSet.find(firstPattern) == _variableLabelSet.end())
               ? firstPattern
               : firstPattern + Constants::ARTIFICIAL_SUFFIX;
  };

  for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
    _constraintMatrix[rowIdx].resize(newVariableCount);
    if (_rowInfos[rowIdx].isObjectiveRow())
      continue;

    const auto newVariableIdx = variableCountAtTheStart + rowIdx - 1;
    _constraintMatrix[rowIdx][newVariableIdx] = 1;
    const auto newArtificialLabelStr = newArtificialLabel(newVariableIdx);
    _variableInfos.push_back(VariableInfo{
        newArtificialLabelStr, VariableType::CONTINUOUS, false, true, true});
    _variableLabelSet.insert(newArtificialLabelStr);
  }
}

#endif // GMISOLVER_LINEARPROGRAM_H
