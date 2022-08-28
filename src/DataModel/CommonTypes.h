#ifndef GMISOLVER_COMMONTYPES_H
#define GMISOLVER_COMMONTYPES_H

#include "src/DataModel/EnumTypes.h"

#include <string>
#include <vector>

struct RowInfo {
  bool isEqualityRow() const { return _type == RowType::EQUALITY; }

  std::string _label;
  RowType _type;
};

struct VariableInfo {
  std::string _label;

  std::string typeStr() const {
    using namespace std::string_literals;

    std::string result;

    if (_type == VariableType::INTEGER)
      result += "INT";

    if (_isSlack)
      result += (!result.empty() ? "|"s : ""s) + "SLACK";

    if (_isArtificial)
      result += (!result.empty() ? "|"s : ""s) + "ARTIFICIAL";

    return result;
  }

  VariableType _type;
  bool _isSlack{false};
  bool _isArtificial{false};
};

template <typename T> using Matrix = std::vector<std::vector<T>>;

#endif // GMISOLVER_COMMONTYPES_H
