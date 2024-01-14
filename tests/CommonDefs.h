#ifndef GMISOLVER_COMMONDEFS_H
#define GMISOLVER_COMMONDEFS_H

#include "src/DataModel/EnumTypes.h"

#include <tuple>

#include <absl/flags/declare.h>
#include <absl/flags/flag.h>

ABSL_DECLARE_FLAG(int32_t, obj_value_logging_frequency);
ABSL_DECLARE_FLAG(int32_t, reinversion_frequency);
ABSL_DECLARE_FLAG(std::vector<SimplexTableauType>, simplex_tableau_types);
ABSL_DECLARE_FLAG(ValidateSimplexOption, validate_simplex_option);
ABSL_DECLARE_FLAG(bool, dual_gomory_remove_only_slack_cuts_with_positive_value);

template <typename... Ts> struct TypeTuple {
  using types = std::tuple<Ts...>;
};

#endif // GMISOLVER_COMMONDEFS_H
