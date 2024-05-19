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
ABSL_FLAG(SlackCutRemovalCondition, slack_cut_removal_condition,
          SlackCutRemovalCondition::ALWAYS,
          "Condition for removing slack cuts");
ABSL_FLAG(bool, extended_statistics, false,
          "Print extended statistics about particular type of optimization");
ABSL_FLAG(int, cut_round_limit, 10,
          "Number of cut rounds during Gomory algorithm");