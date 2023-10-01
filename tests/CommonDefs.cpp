#include "tests/CommonDefs.h"

ABSL_FLAG(
    int32_t, obj_value_logging_frequency, 0,
    "Current objective value should be logged every nth iteration of simplex");
ABSL_FLAG(int32_t, reinversion_frequency, 60,
          "Basis matrix should be reinverted every nth iteration of simplex");
ABSL_FLAG(SimplexTableauType, simplex_tableau_type,
          SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE,
          "Simplex tableau type");
ABSL_FLAG(ValidateSimplexOption, validate_simplex_option,
          ValidateSimplexOption::DONT_VALIDATE,
          "Validate simplex implementations");