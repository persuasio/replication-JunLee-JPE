# ************
# * SCRIPT:   tableD1.R
# * PURPOSE:  Creates Table D1
# *
# * ACKNOWLEDGMENT
# *       The original dataset "FoxNewsFinalDataQJE" is from
# *       DellaVigna and Kaplan (QJE, 2007).
# *
# * NOTES
# *       An alternative measure of the dependent variable is used:
# *       Republican presidential votes as a share of the voting-age population.
# ************

library(haven)
library(dplyr)
library(tibble)
library(fixest)
library(tinytable)
library(here)

output_dir <- here::here("output")
data_path <- here::here(
  "data",
  "DellaVignaKaplan2007",
  "FoxNewsFinalDataQJE.dta"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ",
    data_path,
    call. = FALSE
  )
}

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}


# Original Stata source: trace_tableD1.txt.
# Uses DellaVigna and Kaplan data and fixed-effect regressions to construct
# persuasion-rate bounds. Stata areg is translated to fixest::feols.

weighted_mean_positive <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0

  if (!any(keep)) {
    stop("No observations with nonmissing values and positive weights.", call. = FALSE)
  }

  weighted.mean(x[keep], w[keep])
}

fox <- read_dta(data_path) %>%
  filter(sample12000 == 1) %>%
  mutate(
    new_y2000 = reppresvotes2000 / pop18p2000,
    new_y1996 = reppresvotes1996 / pop18p1996
  )

# Census controls
contrcens2000 <- c(
  "pop2000", "hs2000", "hsp2000", "college2000", "male2000",
  "black2000", "hisp2000", "empl2000", "unempl2000",
  "married2000", "income2000", "urban2000"
)

contrcens00m90 <- c(
  "pop00m90", "hs00m90", "hsp00m90", "college00m90", "male00m90",
  "black00m90", "hisp00m90", "empl00m90", "unempl00m90",
  "married00m90", "income00m90", "urban00m90"
)

# Stata shorthand:
# poptot2000d2-poptot2000d10 noch2000d2-noch2000d10
contrcbl2000 <- c(
  paste0("poptot2000d", 2:10),
  paste0("noch2000d", 2:10)
)

controls <- c(
  contrcens2000,
  contrcens00m90,
  contrcbl2000
)

# Reproduce one column of Stata Table D1.
make_tableD1 <- function(data, fixed_effect) {

  y_formula <- as.formula(
    paste(
      "new_y2000 ~",
      paste(c("new_y1996", "foxnews2000", controls), collapse = " + "),
      "|",
      fixed_effect
    )
  )

  y_fit <- feols(
    y_formula,
    data = data,
    weights = ~ pop18p1996,
    cluster = ~ account2000
  )

  coeff_y <- unname(coef(y_fit)[["foxnews2000"]])
  yhat <- predict(y_fit, newdata = data)

  yhat1 <- clip01(
    yhat + coeff_y - coeff_y * data$foxnews2000
  )

  yhat0 <- clip01(
    yhat - coeff_y * data$foxnews2000
  )

  thetahat_num <- yhat1 - yhat0
  thetahat_den <- pmax(1 - yhat0, 1e-8)

  avg_num <- weighted_mean_positive(
    thetahat_num,
    data$auddiaryScar
  )

  avg_den <- weighted_mean_positive(
    thetahat_den,
    data$auddiaryScar
  )

  apr_lb <- avg_num / avg_den

  t_formula <- as.formula(
    paste(
      "foxanyScar ~",
      paste(c("foxnews2000", controls), collapse = " + "),
      "|",
      fixed_effect
    )
  )

  t_fit <- feols(
    t_formula,
    data = data,
    weights = ~ auddiaryScar,
    cluster = ~ account2000
  )

  coeff_t <- unname(coef(t_fit)[["foxnews2000"]])
  ehat <- predict(t_fit, newdata = data)

  ehat1 <- clip01(
    ehat + coeff_t - coeff_t * data$foxnews2000
  )

  ehat0 <- clip01(
    ehat - coeff_t * data$foxnews2000
  )

  ub_num1 <- yhat1 + 1 - ehat1
  ub_num2 <- yhat0 - ehat0
  ub_num <- pmin(1, ub_num1) - pmax(0, ub_num2)
  ub_den <- pmax(1 - pmax(0, ub_num2), 1e-8)

  apr_ub <- weighted_mean_positive(
    ub_num,
    data$auddiaryScar
  ) / weighted_mean_positive(
    ub_den,
    data$auddiaryScar
  )

  late_weight <- ehat1 - ehat0
  theta_local_den <- pmin(thetahat_den, late_weight)

  local_den <- weighted_mean_positive(
    theta_local_den,
    data$auddiaryScar
  )

  tibble(
    estimand = c(
      "APR (LB)",
      "APR (UB)",
      "LPR (LB)",
      "LPR (UB)"
    ),
    value = c(
      apr_lb,
      apr_ub,
      avg_num / local_den,
      1
    )
  )
}

district_results <- make_tableD1(
  data = fox,
  fixed_effect = "diststate"
)

county_results <- make_tableD1(
  data = fox,
  fixed_effect = "countystate"
)

# Raw numeric matrix in the same row/column orientation as the Stata table.
tableD1_matrix <- district_results %>%
  select(estimand, district = value) %>%
  left_join(
    county_results %>%
      select(estimand, county = value),
    by = "estimand"
  ) %>%
  transmute(
    Estimand = estimand,
    `U.S. House district fixed effects` = district,
    `County fixed effects` = county
  )

# Display version: Stata uses three decimal places and does not multiply by 100.
tableD1_data <- tableD1_matrix %>%
  mutate(
    across(
      -Estimand,
      ~ sprintf("%.3f", .x)
    )
  )

stopifnot(
  nrow(tableD1_matrix) == 4L,
  identical(
    tableD1_matrix$Estimand,
    c("APR (LB)", "APR (UB)", "LPR (LB)", "LPR (UB)")
  )
)

# Save the results using the same output convention as table2.R.
write.csv(
  tableD1_matrix,
  file = file.path(
    output_dir,
    "tableD1_raw.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  tableD1_data,
  file = file.path(
    output_dir,
    "tableD1_display.csv"
  ),
  row.names = FALSE,
  na = ""
)

tableD1_tex <- tt(
  tableD1_data,
  caption = "Table D1. Persuasion Rates: Fox News Effects"
)

print(tableD1_tex)

save_tt(
  tableD1_tex,
  output = file.path(
    output_dir,
    "tableD1.tex"
  ),
  overwrite = TRUE
)
