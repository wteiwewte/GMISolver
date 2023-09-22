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
  else if (string == "UI")
    return BoundType::UPPER_BOUND_INTEGER;

  return std::nullopt;
}

std::set<size_t> allowedLinesCount(const BoundType boundType) {
  switch (boundType) {
  case BoundType::LOWER_BOUND:
    return {4};
  case BoundType::UPPER_BOUND:
    return {4};
  case BoundType::FIXED_VARIABLE:
    return {4};
  case BoundType::BINARY_VARIABLE:
    return {3, 4};
  case BoundType::FREE_VARIABLE:
    return {3};
  case BoundType::LOWER_BOUND_MINUS_INF:
    return {3};
  case BoundType::UPPER_BOUND_INTEGER:
    return {4};
  }
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

std::string
lpOptimizationResultToStr(const LPOptimizationResult lpOptimizationResult) {
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
  case LPOptimizationResult::FAILED_VALIDATION:
    return "FAILED_VALIDATION";
  case LPOptimizationResult::UNKNOWN:
    return "UNKNOWN";
  }

  return "";
}

bool AbslParseFlag(absl::string_view text, ValidateSimplex *validateSimplex,
                   std::string *error) {
  if (text == "yes") {
    *validateSimplex = ValidateSimplex::YES;
    return true;
  }
  if (text == "no") {
    *validateSimplex = ValidateSimplex::NO;
    return true;
  }
  *error = "unknown value for enumeration";
  return false;
}

// AbslUnparseFlag converts from an OutputMode to a string.
// Must be in same namespace as OutputMode.

// Returns a textual flag value corresponding to the OutputMode `mode`.
std::string AbslUnparseFlag(ValidateSimplex validateSimplex) {
  switch (validateSimplex) {
  case ValidateSimplex::YES:
    return "yes";
  case ValidateSimplex::NO:
    return "no";
  default:
    return "unknown";
  }
}