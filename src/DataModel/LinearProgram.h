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
  void convertToStandardForm();
  std::optional<LinearProgram<T>>
  dualProgram(const bool addObjValueFirstVar) const;
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

private:
  friend struct MpsReader<T>;
  template <typename U, typename ComparisonTraitsU> friend class SimplexTableau;
  template <typename U, typename ComparisonTraitsU>
  friend class SimplexValidator;

  std::string basicInformationStr() const;
  void logGeneralInformation() const;
  bool areAllCoefficientsInteger(const int rowIdx) const;

  std::string _name;
  std::vector<RowInfo> _rowInfos;

  std::vector<VariableInfo> _variableInfos;
  boost::dynamic_bitset<> _isVariableFreeBitset;
  std::vector<std::optional<T>> _variableLowerBounds;
  std::vector<std::optional<T>> _variableUpperBounds;
  std::set<std::string> _variableLabelSet;

  std::vector<T> _objective;
  RowInfo _objectiveInfo;
  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;
};

#endif // GMISOLVER_LINEARPROGRAM_H
