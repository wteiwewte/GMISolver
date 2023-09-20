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
  else if (string == "RANGES")
    return SectionType::RANGES;
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
    return BoundType::FIXED_VARIABLE;
  else if (string == "BV")
    return BoundType::BINARY_VARIABLE;
  else if (string == "FR")
    return BoundType::FREE_VARIABLE;
  else if (string == "MI")
    return BoundType::LOWER_BOUND_MINUS_INF;

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
  case RowType::UNKNOWN:
    return "UNKNOWN";
  }

  return "";
}

std::string primalSimplexColumnPivotRuleToStr(
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule) {
  switch (primalSimplexColumnPivotRule) {
  case PrimalSimplexColumnPivotRule::FIRST_ELIGIBLE:
    return "FIRST_ELIGIBLE";
  case PrimalSimplexColumnPivotRule::BIGGEST_ABSOLUTE_REDUCED_COST:
    return "BIGGEST_ABSOLUTE_REDUCED_COST";
  }

  return "";
}

std::string dualSimplexRowPivotRuleToStr(
    const DualSimplexRowPivotRule dualSimplexRowPivotRule) {
  switch (dualSimplexRowPivotRule) {
  case DualSimplexRowPivotRule::FIRST_ELIGIBLE:
    return "FIRST_ELIGIBLE";
  case DualSimplexRowPivotRule::BIGGEST_BOUND_VIOLATION:
    return "BIGGEST_BOUND_VIOLATION";
  }

  return "";
}

std::string lpOptimizationResultToStr(
    const LPOptimizationResult lpOptimizationResult)
{
  switch (lpOptimizationResult) {
  case LPOptimizationResult::UNBOUNDED:
    return "UNBOUNDED";
  case LPOptimizationResult::INFEASIBLE:
    return "INFEASIBLE";
  case LPOptimizationResult::INFEASIBLE_OR_UNBDUNDED:
    return "INFEASIBLE_OR_UNBDUNDED";
  case LPOptimizationResult::BOUNDED_AND_FEASIBLE:
    return "BOUNDED_AND_FEASIBLE";
  case LPOptimizationResult::FAILED_REINVERSION:
    return "FAILED_REINVERSION";
  case LPOptimizationResult::REACHED_ITERATION_LIMIT:
    return "REACHED_ITERATION_LIMIT";
  case LPOptimizationResult::COULD_NOT_LOAD:
    return "COULD_NOT_LOAD";
  case LPOptimizationResult::UNKNOWN:
    return "UNKNOWN";
  }

  return "";
}

