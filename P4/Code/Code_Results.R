# sim_results <- read.csv(#path to sim_results.csv)
library(dplyr)


sim_results_compiled <- sim_results01 |>
  group_by(
    case,method
  ) |>
  summarize(
    mean_bias <- mean(bias),
    coverage_prop <- sum(coverage) / 100
  )


# summarise(
#   n_cases  = n(),
#   mean_age = mean(age_years, na.rm=T),
#   max_age  = max(age_years, na.rm=T),
#   min_age  = min(age_years, na.rm=T),
#   n_males  = sum(gender == "m", na.rm=T))