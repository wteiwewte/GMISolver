#ifndef GMISOLVER_ENUMTYPES_H
#define GMISOLVER_ENUMTYPES_H

#include <cstdint>
#include <optional>
#include <string>

enum class SectionType : int8_t {
  UNDEFINED = 0,
  NAME = 1,
  ROWS,
  COLUMNS,
  RHS,
  BOUNDS,
  END
};

enum class RowType : char {
  EQUALITY = 'E',
  LESS_THAN_OR_EQUAL = 'L',
  GREATER_THAN_OR_EQUAL = 'G',
  OBJECTIVE = 'N'
};

enum class VariableType : int8_t { CONTINUOUS = 0, INTEGER };

enum class BoundType : int8_t {
  LOWER_BOUND = 0,
  UPPER_BOUND,
  FREE_VARIABLE,
  BINARY_VARIABLE
};

enum class LPOptimizationResult : int8_t {
  UNBOUNDED = 0,
  INFEASIBLE,
  BOUNDED_AND_FEASIBLE,
  FAILED_REINVERSION
};

enum class PrimalSimplexColumnPivotRule : int8_t
{
  FIRST_ELIGIBLE = 0,
  BIGGEST_ABSOLUTE_REDUCED_COST
};

enum class DualSimplexRowPivotRule : int8_t
{
  FIRST_ELIGIBLE = 0,
  BIGGEST_BOUND_VIOLATION
};

std::optional<SectionType> stringToSectionType(const std::string &string);
std::optional<RowType> stringToRowType(const std::string &string);
std::optional<BoundType> stringToBoundType(const std::string &string);
std::string rowTypeToStr(const RowType rowType);
std::string primalSimplexColumnPivotRuleToStr(const PrimalSimplexColumnPivotRule primalSimplexColumnPivotRule);
std::string dualSimplexRowPivotRuleToStr(const DualSimplexRowPivotRule dualSimplexRowPivotRule);

#endif // GMISOLVER_ENUMTYPES_H
