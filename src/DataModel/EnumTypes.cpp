#include "src/DataModel/EnumTypes.h"

#include "src/Util/SpdlogHeader.h"

#include <absl/flags/marshalling.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_split.h>

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
  else if (string == "LI")
    return BoundType::LOWER_BOUND_INTEGER;

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
  case BoundType::LOWER_BOUND_INTEGER:
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

std::string lexicographicReoptTypeToStr(
    const LexicographicReoptType lexicographicReoptType) {
  switch (lexicographicReoptType) {
  case LexicographicReoptType::MIN:
    return "MIN";
  case LexicographicReoptType::MAX:
    return "MAX";
  }

  return "";
}

std::string
simplexTableauTypeToStr(const SimplexTableauType simplexTableauType) {
  switch (simplexTableauType) {
  case SimplexTableauType::FULL:
    return "FULL";
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return "REVISED_BASIS_MATRIX_INVERSE";
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE:
    return "REVISED_PRODUCT_FORM_OF_INVERSE";
  }

  return "";
}

std::string matrixRepresentationTypeToStr(
    const MatrixRepresentationType representationType) {
  switch (representationType) {
  case MatrixRepresentationType::NORMAL:
    return "NORMAL";
  case MatrixRepresentationType::SPARSE:
    return "SPARSE";
  }

  return "";
}

bool AbslParseFlag(absl::string_view text,
                   ValidateSimplexOption *validateSimplexOption,
                   std::string *error) {
  if (text == "dont_validate") {
    *validateSimplexOption = ValidateSimplexOption::DONT_VALIDATE;
    return true;
  }
  if (text == "validate_and_stop_on_error") {
    *validateSimplexOption = ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR;
    return true;
  }
  if (text == "validate_and_dont_stop_on_error") {
    *validateSimplexOption =
        ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR;
    return true;
  }
  *error = "unknown value for enumeration";
  return false;
}

std::string AbslUnparseFlag(ValidateSimplexOption validateSimplexOption) {
  switch (validateSimplexOption) {
  case ValidateSimplexOption::DONT_VALIDATE:
    return "dont_validate";
  case ValidateSimplexOption::VALIDATE_AND_STOP_ON_ERROR:
    return "validate_and_stop_on_error";
  case ValidateSimplexOption::VALIDATE_AND_DONT_STOP_ON_ERROR:
    return "validate_and_dont_stop_on_error";
  default:
    return "unknown";
  }
}
bool AbslParseFlag(absl::string_view text,
                   SimplexTableauType *simplexTableauType, std::string *error) {
  if (text == "full") {
    *simplexTableauType = SimplexTableauType::FULL;
    return true;
  }
  if (text == "revised_basis_matrix_inverse") {
    *simplexTableauType = SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE;
    return true;
  }
  if (text == "revised_product_form_of_inverse") {
    *simplexTableauType = SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE;
    return true;
  }
  *error = "unknown value for enumeration";
  return false;
}

std::string AbslUnparseFlag(SimplexTableauType simplexTableauType) {
  switch (simplexTableauType) {
  case SimplexTableauType::FULL:
    return "full";
  case SimplexTableauType::REVISED_BASIS_MATRIX_INVERSE:
    return "revised_basis_matrix_inverse";
  case SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE:
    return "revised_product_form_of_inverse";
  default:
    return "unknown";
  }
}

bool AbslParseFlag(absl::string_view text,
                   std::vector<SimplexTableauType> *simplexTableauTypes,
                   std::string *error) {
  std::vector<absl::string_view> tokens = absl::StrSplit(text, ',');
  simplexTableauTypes->resize(tokens.size());
  for (int tokenIdx = 0; tokenIdx < tokens.size(); ++tokenIdx) {
    SimplexTableauType &currentSimplexTableauType =
        simplexTableauTypes->operator[](tokenIdx);
    if (!absl::ParseFlag(absl::StripLeadingAsciiWhitespace(tokens[tokenIdx]),
                         &currentSimplexTableauType, error))
      return false;
  }

  return true;
}

std::string
AbslUnparseFlag(const std::vector<SimplexTableauType> &simplexTableauTypes) {
  std::vector<std::string> simplexTableauTypesStr(simplexTableauTypes.size());
  for (int typeIdx = 0; typeIdx < simplexTableauTypes.size(); ++typeIdx)
    simplexTableauTypesStr[typeIdx] =
        AbslUnparseFlag(simplexTableauTypes[typeIdx]);

  return absl::StrJoin(simplexTableauTypesStr, ", ");
}