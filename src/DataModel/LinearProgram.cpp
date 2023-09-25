#include "src/DataModel/LinearProgram.h"

#include "src/Util/LPPrinter.h"

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
         std::all_of(
             _variableLowerBounds.begin(), _variableLowerBounds.end(),
             [](const std::optional<T> &lb) { return lb.has_value(); }) &&
         (_variableUpperBounds.size() == _variableInfos.size()) &&
         std::all_of(_variableUpperBounds.begin(), _variableUpperBounds.end(),
                     [](const std::optional<T> &up) { return up.has_value(); });
}

template <typename T> bool LinearProgram<T>::isPureIP() const {
  return std::all_of(_variableInfos.begin(), _variableInfos.end(),
                     [](const VariableInfo &varInfo) {
                       return varInfo._type == VariableType::INTEGER;
                     });
}

template class LinearProgram<double>;
template class LinearProgram<long double>;