#include "tests/CommonDefs.h"

ABSL_FLAG(
    int32_t, obj_value_logging_frequency, 0,
    "Current objective value should be logged every nth iteration of simplex");
ABSL_FLAG(int32_t, reinversion_frequency, 60,
          "Basis matrix should be reinverted every nth iteration of simplex");
ABSL_FLAG(bool, use_product_form_of_inverse, true,
          "Basis matrix inverse is represented via product form of inverse");
ABSL_FLAG(ValidateSimplexOption, validate_simplex_option,
          ValidateSimplexOption::DONT_VALIDATE,
          "Validate simplex implementations");