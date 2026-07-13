# ************
# * SCRIPT:   tableH1.R
# * PURPOSE:  Creates Table H1
# *
# * ACKNOWLEDGMENT
# *       The original dataset "NTV_Individual_Data.dta" is from
# *       Enikolopov, Petrova, and Zhuravskaya (AER, 2011).
# ************

library(haven)
library(dplyr)
library(tibble)
library(survey)
library(tinytable)
library(here)

output_dir <- here::here("output")
data_path <- here::here(
  "data",
  "EnikolopovPetrovaZhuravskaya2011",
  "NTV_Individual_Data.dta"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ",
    data_path,
    call. = FALSE
  )
}

basic <- c("logpop98", "wage98_ln")
sociodem <- c("male", "age", "educ1", "married", "consump")
controls <- c(sociodem, basic)
parties <- c("Unity", "OVR")

# Reproduce Stata's pweight percentile sufficiently closely for this binary split
weighted_quantile <- function(x, w, probability = 0.5) {
  keep <- is.finite(x) & is.finite(w) & w > 0

  if (!any(keep)) {
    stop("No observations with nonmissing values and positive weights.", call. = FALSE)
  }

  x <- x[keep]
  w <- w[keep]
  ordering <- order(x)
  x <- x[ordering]
  w <- w[ordering]

  cumulative_weight <- cumsum(w) / sum(w)
  x[which(cumulative_weight >= probability)[1L]]
}

make_formula <- function(outcome) {
  reformulate(
    termlabels = controls,
    response = outcome
  )
}

fit_weighted_probit <- function(formula, data) {
  suppressWarnings(
    glm(
      formula = formula,
      data = data,
      weights = kishweig,
      family = quasibinomial(link = "probit"),
      na.action = na.exclude
    )
  )
}

predict_response <- function(model, data) {
  as.numeric(
    predict(
      model,
      newdata = data,
      type = "response"
    )
  )
}

# Stata's svyset [pweight = kishweig] treats observations as PSUs.
svy_weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0

  design <- svydesign(
    ids = ~1,
    weights = ~weight,
    data = data.frame(
      value = x[keep],
      weight = w[keep]
    )
  )

  unname(coef(svymean(~value, design))[[1L]])
}

svy_weighted_mean_se <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0

  design <- svydesign(
    ids = ~1,
    weights = ~weight,
    data = data.frame(
      value = x[keep],
      weight = w[keep]
    )
  )

  unname(SE(svymean(~value, design))[[1L]])
}

ntv <- read_dta(data_path) %>%
  mutate(
    # Match Stata so that a missing vote_Unity or vote_OVR remains 0
    y_vote_Unity = as.integer(
      !is.na(.data$vote_Unity) & .data$vote_Unity == 0
    ),

    # Match Stata so that a missing vote_Unity or vote_OVR remains 0
    y_vote_OVR = as.integer(
      !is.na(.data$vote_OVR) & .data$vote_OVR == 1
    )
  )

required_variables <- c(
  "kishweig",
  "tvmaxtveloss5050powerA",
  "Watches_NTV_1999",
  "vote_Unity",
  "vote_OVR",
  controls
)

missing_variables <- setdiff(required_variables, names(ntv))

if (length(missing_variables) > 0L) {
  stop(
    "Required variables are missing: ",
    paste(missing_variables, collapse = ", "),
    call. = FALSE
  )
}

# Create the binary instrument using the weighted median of signal strength.
IV_p50 <- weighted_quantile(
  x = ntv$tvmaxtveloss5050powerA,
  w = ntv$kishweig,
  probability = 0.5
)

ntv <- ntv %>%
  mutate(
    IV_NTV = as.integer(.data$tvmaxtveloss5050powerA > IV_p50)
  )

# Predict P(Z = 1 | X).
z_fit <- fit_weighted_probit(
  make_formula("IV_NTV"),
  data = ntv
)

ntv$phat_Z1 <- predict_response(z_fit, ntv)

# These treatment regressions appear in the Stata trace even though their
# fitted values are not in the final Table H1 calculations.
t_fit_z1 <- fit_weighted_probit(
  make_formula("Watches_NTV_1999"),
  data = ntv %>% filter(.data$IV_NTV == 1)
)

t_fit_z0 <- fit_weighted_probit(
  make_formula("Watches_NTV_1999"),
  data = ntv %>% filter(.data$IV_NTV == 0)
)

ntv$phat_T_Z1 <- predict_response(t_fit_z1, ntv)
ntv$phat_T_Z0 <- predict_response(t_fit_z0, ntv)

estimate_party <- function(data, party) {
  outcome <- paste0("y_vote_", party)

  fit_z1 <- fit_weighted_probit(
    make_formula(outcome),
    data = data %>% filter(.data$IV_NTV == 1)
  )

  fit_z0 <- fit_weighted_probit(
    make_formula(outcome),
    data = data %>% filter(.data$IV_NTV == 0)
  )

  phat_y_z1 <- predict_response(fit_z1, data)
  phat_y_z0 <- predict_response(fit_z0, data)
  observed_y <- data[[outcome]]

  # Augmented inverse-probability terms, following Stata convention
  w1 <- data$IV_NTV / data$phat_Z1
  w0 <- (1 - data$IV_NTV) / (1 - data$phat_Z1)

  t1 <- w1 * (observed_y - phat_y_z1)
  t0 <- w0 * (observed_y - phat_y_z0)

  beta1 <- svy_weighted_mean(phat_y_z1, data$kishweig)
  beta0 <- svy_weighted_mean(phat_y_z0, data$kishweig)

  lower_bound <- (beta1 - beta0) / (1 - beta0)

  infl_t1 <- t1 + (phat_y_z1 - beta1)
  infl_t0 <- t0 + (phat_y_z0 - beta0)

  influence <- (1 / (1 - beta0)) *
    (
      infl_t1 -
        ((1 - beta1) / (1 - beta0)) * infl_t0
    )

  standard_error <- svy_weighted_mean_se(
    influence,
    data$kishweig
  )

  # A 1.645 critical value gives a one-sided 95% lower confidence limit.
  confidence_lower <- lower_bound - 1.645 * standard_error

  c(
    point_estimate = lower_bound,
    standard_error = standard_error,
    confidence_lower = confidence_lower
  )
}

unity_results <- estimate_party(ntv, "Unity")
ovr_results <- estimate_party(ntv, "OVR")

# Row labels used in both raw and display tables
row_labels <- c(
  "Point estimate of the lower bound",
  "Standard error of the lower bound",
  "One-sided 95% confidence interval for θ_avg"
)

# Raw numeric results
tableH1_raw <- tibble(
  Statistic = row_labels,
  `Not voting for Unity` = c(
    unity_results["point_estimate"],
    unity_results["standard_error"],
    unity_results["confidence_lower"]
  ),
  `Voting for OVR` = c(
    ovr_results["point_estimate"],
    ovr_results["standard_error"],
    ovr_results["confidence_lower"]
  )
)

# Formatted results for display
tableH1_display <- tableH1_raw %>%
  mutate(
    `Not voting for Unity` = c(
      sprintf("%.3f", .data$`Not voting for Unity`[1]),
      sprintf("%.3f", .data$`Not voting for Unity`[2]),
      sprintf("[%.3f,1]", .data$`Not voting for Unity`[3])
    ),
    `Voting for OVR` = c(
      sprintf("%.3f", .data$`Voting for OVR`[1]),
      sprintf("%.3f", .data$`Voting for OVR`[2]),
      sprintf("[%.3f,1]", .data$`Voting for OVR`[3])
    )
  )

# Save raw and formatted CSV files
write.csv(
  tableH1_raw,
  file = file.path(output_dir, "tableH1_raw.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  tableH1_display,
  file = file.path(output_dir, "tableH1_display.csv"),
  row.names = FALSE,
  na = ""
)

# Use LaTeX math only in the LaTeX version
tableH1_latex <- tableH1_display
tableH1_latex$Statistic[3] <-
  "One-sided 95\\% confidence interval for $\\theta_{\\mathrm{avg}}$"

tableH1_tex <- tt(
  tableH1_latex,
  caption = "TABLE H1. Persuasion Rates: NTV Effects Using a Binary Instrument"
) |>
  style_tt(
    j = 1,
    align = "l"
  ) |>
  style_tt(
    j = 2:3,
    align = "c"
  ) |>
  group_tt(
    j = list(
      "(1)" = 2,
      "(2)" = 3
    )
  )

print(tableH1_tex)

save_tt(
  tableH1_tex,
  output = file.path(output_dir, "tableH1.tex"),
  overwrite = TRUE
)
