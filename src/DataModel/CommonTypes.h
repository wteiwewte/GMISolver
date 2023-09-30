#ifndef GMISOLVER_COMMONTYPES_H
#define GMISOLVER_COMMONTYPES_H

#include "src/DataModel/EnumTypes.h"

#include <string>
#include <vector>

struct RowInfo {
  bool isEqualityRow() const { return _type == RowType::EQUALITY; }

  std::string _label;
  RowType _type = RowType::UNKNOWN;
};

struct VariableInfo {
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

  std::string _label;
  VariableType _type;
  bool _isSlack{false};
  bool _isArtificial{false};
  bool _isFixed{false};
  bool _isObjectiveVar{false};
  std::optional<int> _cutRowIdx;
};

template <typename T> struct PivotData {
  const int _leavingRowIdx;
  const int _enteringColumnIdx;
  const T _pivotingTermInverse;
};

template <typename T> struct PivotRowData {
  std::optional<int> _pivotRowIdx;
  std::optional<T> _minRatio;
  bool _departingIdxBecomesLowerBound{false};
  bool _noBasisChangeNeeded{false};
};

#endif // GMISOLVER_COMMONTYPES_H
