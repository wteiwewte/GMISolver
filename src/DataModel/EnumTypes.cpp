#include "src/DataModel/EnumTypes.h"

std::optional<SectionType> stringToSectionType(const std::string &string) {
  if (string == "NAME")
    return SectionType::NAME;
  else if (string == "ROWS")
    return SectionType::ROWS;
  else if (string == "COLUMNS")
    return SectionType::COLUMNS;
  else if (string == "RHS")
    return SectionType::RHS;
  else if (string == "BOUNDS")
    return SectionType::BOUNDS;
  else if (string == "ENDATA")
    return SectionType::END;

  return std::nullopt;
}
std::optional<RowType> stringToRowType(const std::string &string) {
  if (string[0] == static_cast<char>(RowType::EQUALITY))
    return RowType::EQUALITY;
  if (string[0] == static_cast<char>(RowType::LESS_THAN_OR_EQUAL))
    return RowType::LESS_THAN_OR_EQUAL;
  if (string[0] == static_cast<char>(RowType::GREATER_THAN_OR_EQUAL))
    return RowType::GREATER_THAN_OR_EQUAL;
  if (string[0] == static_cast<char>(RowType::OBJECTIVE))
    return RowType::OBJECTIVE;

  return std::nullopt;
}

std::optional<BoundType> stringToBoundType(const std::string &string) {
  if (string == "LO")
    return BoundType::LOWER_BOUND;
  else if (string == "UP")
    return BoundType::UPPER_BOUND;
  else if (string == "FX")
      return BoundType::FREE_VARIABLE;
  else if (string == "BV")
    return BoundType::BINARY_VARIABLE;

  return std::nullopt;
}

std::string rowTypeToStr(const RowType rowType) {
  switch (rowType) {
  case RowType::EQUALITY:
    return "=";
  case RowType::LESS_THAN_OR_EQUAL:
    return "<=";
  case RowType::GREATER_THAN_OR_EQUAL:
    return ">=";
  case RowType::OBJECTIVE:
    return "OBJ";
  }

  return "";
}

std::string primalSimplexColumnPivotRuleToStr(const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule)
{
  switch (primalSimplexColumnPivotRule) {
  case PrimalSimplexColumnPivotRule::FIRST_ELIGIBLE:
    return "FIRST_ELIGIBLE";
  case PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST:
    return "BIGGEST_ABSOLUTE_REDUCED_COST";
  }

  return "";
}

std::string dualSimplexRowPivotRuleToStr(const DualSimplexRowPivotRule dualSimplexRowPivotRule)
{
  switch (dualSimplexRowPivotRule) {
  case DualSimplexRowPivotRule::FIRST_ELIGIBLE:
    return "FIRST_ELIGIBLE";
  case DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION:
    return "BIGGEST_BOUND_VIOLATION";
  }

  return "";
}
