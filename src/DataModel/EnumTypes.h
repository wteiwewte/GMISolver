#ifndef GMISOLVER_ENUMTYPES_H
#define GMISOLVER_ENUMTYPES_H

#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <variant>

#include <absl/strings/string_view.h>

enum class SectionType : int8_t {
  UNDEFINED = 0,
  NAME,
  ROWS,
  COLUMNS,
  RHS,
  RANGES,
  BOUNDS,
  END
};

enum class RowType : char {
  EQUALITY = 'E',
  LESS_THAN_OR_EQUAL = 'L',
  GREATER_THAN_OR_EQUAL = 'G',
  OBJECTIVE = 'N',
  UNKNOWN = 'U'
};

enum class VariableType : int8_t { CONTINUOUS = 0, INTEGER };

enum class BoundType : int8_t {
  LOWER_BOUND = 0,
  UPPER_BOUND,
  FIXED_VARIABLE,
  BINARY_VARIABLE,
  FREE_VARIABLE,
  LOWER_BOUND_MINUS_INF,
  UPPER_BOUND_INTEGER,
  LOWER_BOUND_INTEGER,
  UPPER_BOUND_PLUS_INF
};

enum class LPOptimizationResult : int8_t {
  UNBOUNDED = 0,
  INFEASIBLE,
  INFEASIBLE_OR_UNBDUNDED,
  BOUNDED_AND_FEASIBLE,
  FAILED_REINVERSION,
  REACHED_ITERATION_LIMIT,
  COULD_NOT_LOAD,
  FAILED_VALIDATION,
  TOO_SMALL_PROGRESS,
  UNKNOWN
};

enum class PrimalSimplexColumnPivotRule : int8_t {
  FIRST_ELIGIBLE = 0,
  BIGGEST_ABSOLUTE_REDUCED_COST
};

enum class DualSimplexRowPivotRule : int8_t {
  FIRST_ELIGIBLE = 0,
  BIGGEST_BOUND_VIOLATION
};

enum class LPOptimizationType : int8_t {
  LINEAR_RELAXATION = 0,
  INTEGER_PROGRAM,
  MIXED_INTEGER_PROGRAM
};

enum class SimplexType : int8_t {
  PRIMAL = 0,
  DUAL,
};

enum class SimplexTableauType : int8_t {
  FULL = 0,
  REVISED_BASIS_MATRIX_INVERSE,
  REVISED_PRODUCT_FORM_OF_INVERSE
};

enum class MatrixRepresentationType : int8_t {
  NORMAL = 0,
  SPARSE,
};

enum class LexicographicReoptType : int8_t {
  MIN = 0,
  MAX,
};

enum class GomoryCutChoosingRule : int8_t { FIRST = 0, ALL, RANDOM };

enum class ValidateSimplexOption : int8_t {
  DONT_VALIDATE = 0,
  VALIDATE_AND_STOP_ON_ERROR,
  VALIDATE_AND_DONT_STOP_ON_ERROR,
  UNKNOWN
};

enum class SlackCutRemovalCondition : int8_t {
  ALWAYS = 0,
  ONLY_WHEN_SLACK_VAR_IS_POSITIVE
};

enum class SaveLexSolution : int8_t { YES = 0, NO };

enum class PrimalPhase : int8_t { ONE = 0, TWO };
enum class DualPhase : int8_t { ONE = 0, TWO };
using SimplexPhase = std::variant<PrimalPhase, DualPhase>;

enum class PrintSimplexOptSummary : int8_t { YES = 0, NO };
enum class AddObjectiveRelatedVar : int8_t { YES = 0, NO };
enum class IsDualProgramOptimized : int8_t { YES = 0, NO };
enum class DoVarsMustBeNonnegative : int8_t { YES = 0, NO };
enum class AllBoundsMustBeSpecified : int8_t { YES = 0, NO };
enum class IsPrimalCuttingPlanes : int8_t { YES = 0, NO };

bool AbslParseFlag(absl::string_view text,
                   ValidateSimplexOption *validateSimplexOption,
                   std::string *error);
std::string AbslUnparseFlag(ValidateSimplexOption validateSimplexOption);
bool AbslParseFlag(absl::string_view text,
                   SimplexTableauType *simplexTableauType, std::string *error);
std::string AbslUnparseFlag(SimplexTableauType simplexTableauType);
bool AbslParseFlag(absl::string_view text,
                   SlackCutRemovalCondition *slackCutRemovalCondition,
                   std::string *error);
std::string AbslUnparseFlag(SlackCutRemovalCondition slackCutRemovalCondition);
bool AbslParseFlag(absl::string_view text,
                   std::vector<SimplexTableauType> *simplexTableauTypes,
                   std::string *error);
std::string
AbslUnparseFlag(const std::vector<SimplexTableauType> &simplexTableauTypes);

std::optional<SectionType> stringToSectionType(const std::string &string);
std::optional<RowType> stringToRowType(const std::string &string);
std::optional<BoundType> stringToBoundType(const std::string &string);
std::set<size_t> allowedLinesCount(const BoundType boundType);
std::string rowTypeToStr(const RowType rowType);
std::string primalSimplexColumnPivotRuleToStr(
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule);
std::string dualSimplexRowPivotRuleToStr(
    const DualSimplexRowPivotRule dualSimplexRowPivotRule);
std::string
lpOptimizationResultToStr(const LPOptimizationResult lpOptimizationResult);
std::string lexicographicReoptTypeToStr(
    const LexicographicReoptType lexicographicReoptType);
std::string
simplexTableauTypeToStr(const SimplexTableauType simplexTableauType);
std::string matrixRepresentationTypeToStr(
    const MatrixRepresentationType representationType);

#endif // GMISOLVER_ENUMTYPES_H
