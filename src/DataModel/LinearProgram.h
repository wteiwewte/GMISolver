#ifndef GMISOLVER_LINEARPROGRAM_H
#define GMISOLVER_LINEARPROGRAM_H

#include "src/DataModel/CommonConstants.h"
#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/EnumTypes.h"
#include "src/DataModel/MatrixTypes.h"

#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

template <typename T> struct MpsReader;

template <typename T> class LinearProgram {
public:
  void convertAllConstraintsToEquality();
  LinearProgram<T> convertToStandardForm() const;
  std::optional<LinearProgram<T>>
  dualProgram(const AddObjectiveRelatedVar addObjectiveRelatedVar) const;
  size_t getOriginalVariablesCount() const;

  std::string toString() const;
  std::string toStringLpSolveFormat() const;
  void logGeneralInfo() const;

  const std::vector<RowInfo> &getRowInfos() const { return _rowInfos; }
  const std::vector<VariableInfo> &getVariableInfos() const {
    return _variableInfos;
  }
  const std::vector<T> &getObjective() const { return _objective; }
  const std::string &getName() const { return _name; }

  bool checkIfAllBoundsAreSpeficied() const;
  bool isPureIP(const bool checkObjVar = true) const;
  bool allCoefficientsAreIntegers() const;
  bool allVariablesAreNonnegative() const;
  T objectiveConstant() const { return _objectiveConstant; }

private:
  friend struct MpsReader<T>;
  template <typename U, typename ComparisonTraitsU> friend class SimplexTableau;
  template <typename U, typename ComparisonTraitsU>
  friend class SimplexValidator;

  std::string basicInformationStr() const;
  void logGeneralInformation() const;
  bool areAllCoefficientsInteger(const int rowIdx) const;

  static void
  addVariablesFromOriginLP(const LinearProgram<T> &originLPWithEquality,
                           LinearProgram<T> &standardFormLP);
  static void negateVariablesIfNeeded(LinearProgram<T> &standardFormLP);
  static void makeAllVariablesNonnegative(LinearProgram<T> &standardFormLP);

  std::string _name;
  std::vector<RowInfo> _rowInfos;

  std::vector<VariableInfo> _variableInfos;
  std::vector<std::optional<T>> _variableLowerBounds;
  std::vector<std::optional<T>> _variableUpperBounds;
  std::set<std::string> _variableLabelSet;

  std::vector<T> _objective;
  RowInfo _objectiveInfo;
  T _objectiveConstant{0};
  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;
};

#endif // GMISOLVER_LINEARPROGRAM_H
