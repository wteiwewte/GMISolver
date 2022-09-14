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

template class LinearProgram<double>;