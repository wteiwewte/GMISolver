#ifndef GMISOLVER_COMMONDEFS_H
#define GMISOLVER_COMMONDEFS_H

#include "src/DataModel/EnumTypes.h"

#include <tuple>

#include <absl/flags/declare.h>
#include <absl/flags/flag.h>

ABSL_DECLARE_FLAG(int32_t, obj_value_logging_frequency);
ABSL_DECLARE_FLAG(int32_t, reinversion_frequency);
ABSL_DECLARE_FLAG(bool, use_product_form_of_inverse);
ABSL_DECLARE_FLAG(ValidateSimplex, validate_simplex);

template <typename... Ts> struct TypeTuple {
  using types = std::tuple<Ts...>;
};

#endif // GMISOLVER_COMMONDEFS_H
