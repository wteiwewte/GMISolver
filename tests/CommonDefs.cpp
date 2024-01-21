#include "tests/CommonDefs.h"

ABSL_FLAG(
    int32_t, obj_value_logging_frequency, 0,
    "Current objective value should be logged every nth iteration of simplex");
ABSL_FLAG(int32_t, reinversion_frequency, 60,
          "Basis matrix should be reinverted every nth iteration of simplex");
ABSL_FLAG(std::vector<SimplexTableauType>, simplex_tableau_types,
          std::vector<SimplexTableauType>(
              {SimplexTableauType::REVISED_PRODUCT_FORM_OF_INVERSE}),
          "Simplex tableau types");
ABSL_FLAG(ValidateSimplexOption, validate_simplex_option,
          ValidateSimplexOption::DONT_VALIDATE,
          "Validate simplex implementations");
ABSL_FLAG(bool, dual_gomory_remove_only_slack_cuts_with_positive_value, false,
          "Dual Gomory - remove slack cuts only when its corresponding "
          "variable has value greater than 0");
ABSL_FLAG(bool, extended_statistics, false,
          "Print extended statistics about particular type of optimization");