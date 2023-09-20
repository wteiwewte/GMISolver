#ifndef GMISOLVER_ENUMTYPES_H
#define GMISOLVER_ENUMTYPES_H

#include <cstdint>
#include <optional>
#include <string>

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
  LOWER_BOUND_MINUS_INF
};

enum class LPOptimizationResult : int8_t {
  UNBOUNDED = 0,
  INFEASIBLE,
  INFEASIBLE_OR_UNBDUNDED,
  BOUNDED_AND_FEASIBLE,
  FAILED_REINVERSION,
  REACHED_ITERATION_LIMIT,
  COULD_NOT_LOAD,
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

std::optional<SectionType> stringToSectionType(const std::string &string);
std::optional<RowType> stringToRowType(const std::string &string);
std::optional<BoundType> stringToBoundType(const std::string &string);
std::string rowTypeToStr(const RowType rowType);
std::string primalSimplexColumnPivotRuleToStr(
    const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule);
std::string dualSimplexRowPivotRuleToStr(
    const DualSimplexRowPivotRule dualSimplexRowPivotRule);
std::string lpOptimizationResultToStr(
    const LPOptimizationResult lpOptimizationResult);

#endif // GMISOLVER_ENUMTYPES_H
