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

template <typename T> struct MpsReader;

template <typename T> class LinearProgram {
public:
  std::string toString() const;
  std::string toStringLpSolveFormat() const;

  const std::vector<RowInfo> &getRowInfos() const { return _rowInfos; }
  const std::vector<VariableInfo> &getVariableInfos() const {
    return _variableInfos;
  }
  const Matrix<T> &getConstraintMatrix() const { return _constraintMatrix; }
  const std::vector<T> &getRightHandSides() const { return _rightHandSides; }
  const std::vector<T> &getObjective() const { return _objective; }

private:
  friend struct MpsReader<T>;
  template <typename U, typename ComparisonTraitsU>
  friend class SimplexTableau;

  std::string basicInformationStr() const;

  std::string _name;
  std::vector<RowInfo> _rowInfos;

  std::vector<VariableInfo> _variableInfos;
  std::vector<std::optional<T>> _variableLowerBounds;
  std::vector<std::optional<T>> _variableUpperBounds;
  std::set<std::string> _variableLabelSet;

  std::vector<T> _objective;
  RowInfo _objectiveInfo;
  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;
};

#endif // GMISOLVER_LINEARPROGRAM_H
