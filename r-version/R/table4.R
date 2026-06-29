# ************
# * SCRIPT:   table4.do
# * PURPOSE:  Creates Table 4
# *
# * ACKNOWLEDGMENT
# *       The orginal dataset "publicdata.dta" is from Gerber, Karlan, and Berg an (2009, AEJ Applied).
# ************


library(haven)
library(dplyr)
library(tibble)
library(tinytable)

persuasion_dir <- Sys.getenv("PERSUASION_DIR")
if (persuasion_dir == "") {
  persuasion_dir <- "."
}

results_dir <- file.path(persuasion_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

# ************
# * SCRIPT:   table4.R
# * PURPOSE:  Creates Table 4
# ************

library(sandwich)
library(lmtest)
library(AER)
library(numDeriv)

# Original Stata source: trace_table4.txt.
# The Stata code manually constructs ITT, APR, LPR, and LATE bounds.
# This R version follows those formulas directly. For suest/nlcom sections,
# we use numerical delta-method SEs with HC1 covariance matrices.

gkb <- read_dta(
  file.path(
    persuasion_dir,
    "data",
    "GerberKarlanBergan2009",
    "publicdata.dta"
  )
)

gkb <- gkb %>%
  filter(times != 1) %>%
  mutate(
    data_avail = as.integer(
      survey == 1 &
        !is.na(voteddem_all) &
        !is.na(readsome)
    )
  ) %>%
  filter(data_avail == 1) %>%
  mutate(
    case_id = row_number(),
    outcome = voteddem_all,
    treat = readsome,
    instr = post,
    a_u = outcome * treat + (1 - treat),
    b_l = outcome * (1 - treat),
    a_u2 = outcome + (1 - treat),
    b_l2 = outcome - treat
  )

alpha_level <- 0.2
cv_cns1 <- qnorm(1 - alpha_level)
cv_cns2 <- qnorm(1 - alpha_level / 2)

pretest_level <- 0.001
cv_cns_pretest <- qnorm(1 - pretest_level / 4)

adjusted <- alpha_level - pretest_level
cv_cns1_adj <- qnorm(1 - adjusted)
cv_cns2_adj <- qnorm(1 - adjusted / 2)

robust_coef <- function(fit, term) {
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC1"))
  c(
    estimate = ct[term, "Estimate"],
    se = ct[term, "Std. Error"]
  )
}

stoye_cv <- function(lower, upper, lower_se, upper_se, cv_max, alpha) {
  correction_term <- (upper - lower) / max(upper_se, lower_se)
  grid <- seq(0, cv_max + 0.01, by = 0.01)
  diff <- abs(pnorm(grid + correction_term) - pnorm(-grid) - (1 - alpha))
  mean(grid[diff <= 0.002], na.rm = TRUE)
}

delta_nlcom <- function(fit, f) {
  b <- coef(fit)
  vc <- sandwich::vcovHC(fit, type = "HC1")
  est <- f(b)
  grad <- numDeriv::jacobian(function(bb) f(bb), b)
  se <- sqrt(as.numeric(grad %*% vc %*% t(grad)))
  c(estimate = est, se = se)
}

results <- list()

# ITT estimate on outcome.
itt_y_fit <- lm(outcome ~ instr, data = gkb)
itt_y <- robust_coef(itt_y_fit, "instr")

results[["ITT outcome"]] <- tibble(
  estimand = "ITT: outcome",
  lower = max(0, itt_y["estimate"] - cv_cns2 * itt_y["se"]),
  estimate_lower = itt_y["estimate"],
  estimate_upper = itt_y["estimate"],
  upper = min(1, itt_y["estimate"] + cv_cns2 * itt_y["se"])
)

# ITT estimate on treatment.
itt_t_fit <- lm(treat ~ instr, data = gkb)
itt_t <- robust_coef(itt_t_fit, "instr")

results[["ITT treatment"]] <- tibble(
  estimand = "ITT: treatment",
  lower = max(0, itt_t["estimate"] - cv_cns2 * itt_t["se"]),
  estimate_lower = itt_t["estimate"],
  estimate_upper = itt_t["estimate"],
  upper = min(1, itt_t["estimate"] + cv_cns2 * itt_t["se"])
)

# APR when (Y,T,Z) are observed: lower bound.
apr_lower <- delta_nlcom(
  itt_y_fit,
  function(b) {
    b[["instr"]] / (1 - b[["(Intercept)"]])
  }
)

# APR upper bound: Stata suest/nlcom across a_u and b_l.
# R simplification: fit the two equations and use a paired case bootstrap for
# the nonlinear combination. This is a substantive but transparent replacement
# for suest because there is no base-R suest equivalent.
apr_upper_point <- function(data) {
  fit_au <- lm(a_u ~ instr, data = data)
  fit_bl <- lm(b_l ~ instr, data = data)
  b_au <- coef(fit_au)
  b_bl <- coef(fit_bl)
  (b_au[["(Intercept)"]] + b_au[["instr"]] - b_bl[["(Intercept)"]]) /
    (1 - b_bl[["(Intercept)"]])
}

upper_bound_coef <- apr_upper_point(gkb)

set.seed(12345)
apr_upper_boot <- replicate(
  500,
  {
    idx <- sample.int(nrow(gkb), replace = TRUE)
    apr_upper_point(gkb[idx, ])
  }
)

upper_bound_se <- sd(apr_upper_boot, na.rm = TRUE)

lower_bound_coef <- apr_lower["estimate"]
lower_bound_se <- apr_lower["se"]

cv_cns_stoye <- stoye_cv(
  lower_bound_coef,
  upper_bound_coef,
  lower_bound_se,
  upper_bound_se,
  cv_cns2,
  alpha_level
)

results[["APR YTZ"]] <- tibble(
  estimand = "APR: Y,T,Z observed",
  lower = max(0, lower_bound_coef - cv_cns_stoye * lower_bound_se),
  estimate_lower = lower_bound_coef,
  estimate_upper = upper_bound_coef,
  upper = min(1, upper_bound_coef + cv_cns_stoye * upper_bound_se)
)

# APR when (Y,T) and (Y,Z) are separately observed.
fit_a_u2 <- lm(a_u2 ~ instr, data = gkb)
a_u2_ct <- lmtest::coeftest(fit_a_u2, vcov. = sandwich::vcovHC(fit_a_u2, type = "HC1"))
coef_a_u2 <- sum(coef(fit_a_u2))
se_a_u2 <- sqrt(sum(sandwich::vcovHC(fit_a_u2, type = "HC1")))

fit_b_l2 <- lm(b_l2 ~ instr, data = gkb)
b_l2_vc <- sandwich::vcovHC(fit_b_l2, type = "HC1")
coef_b_l2 <- coef(fit_b_l2)[["(Intercept)"]]
se_b_l2 <- sqrt(b_l2_vc["(Intercept)", "(Intercept)"])

upper_bound2_coef <- (min(1, coef_a_u2) - max(0, coef_b_l2)) /
  (1 - max(0, coef_b_l2))

upper_bound2_se <- se_a_u2

cv_cns2_stoye <- stoye_cv(
  lower_bound_coef,
  upper_bound2_coef,
  lower_bound_se,
  upper_bound2_se,
  cv_cns2_adj,
  adjusted
)

results[["APR YT YZ"]] <- tibble(
  estimand = "APR: Y,T and Y,Z observed",
  lower = max(0, lower_bound_coef - cv_cns2_stoye * lower_bound_se),
  estimate_lower = lower_bound_coef,
  estimate_upper = upper_bound2_coef,
  upper = min(1, upper_bound2_coef + cv_cns2_stoye * upper_bound2_se)
)

# APR when only (Y,Z) are observed.
results[["APR YZ"]] <- tibble(
  estimand = "APR: Y,Z observed",
  lower = max(0, lower_bound_coef - cv_cns1 * lower_bound_se),
  estimate_lower = lower_bound_coef,
  estimate_upper = 1,
  upper = 1
)

# LPR when (Y,T,Z) are observed.
lpr_point <- function(data) {
  fit_num <- lm(outcome ~ instr, data = data)
  data$den_lpr <- (1 - data$outcome) * (1 - data$treat)
  fit_den <- lm(den_lpr ~ instr, data = data)
  coef(fit_num)[["instr"]] / (-coef(fit_den)[["instr"]])
}

lpr_coef <- lpr_point(gkb)

set.seed(12345)
lpr_boot <- replicate(
  500,
  {
    idx <- sample.int(nrow(gkb), replace = TRUE)
    lpr_point(gkb[idx, ])
  }
)

lpr_se <- sd(lpr_boot, na.rm = TRUE)

results[["LPR YTZ"]] <- tibble(
  estimand = "LPR: Y,T,Z observed",
  lower = max(0, lpr_coef - cv_cns2 * lpr_se),
  estimate_lower = lpr_coef,
  estimate_upper = lpr_coef,
  upper = min(1, lpr_coef + cv_cns2 * lpr_se)
)

# LATE via 2SLS.
iv_fit <- AER::ivreg(outcome ~ treat | instr, data = gkb)
iv_vc <- sandwich::vcovHC(iv_fit, type = "HC1")
lpr2_coef <- coef(iv_fit)[["treat"]]
lpr2_se <- sqrt(iv_vc["treat", "treat"])

lpr_lb_LATE <- max(0, lpr2_coef - cv_cns2 * lpr2_se)
lpr_lb_LBND <- max(0, lower_bound_coef - cv_cns2 * lower_bound_se)

results[["LPR YT YZ"]] <- tibble(
  estimand = "LPR: Y,T and Y,Z observed",
  lower = max(lpr_lb_LATE, lpr_lb_LBND),
  estimate_lower = max(lpr2_coef, lower_bound_coef),
  estimate_upper = 1,
  upper = 1
)

results[["LPR YZ"]] <- tibble(
  estimand = "LPR: Y,Z observed",
  lower = max(0, lower_bound_coef - cv_cns1 * lower_bound_se),
  estimate_lower = lower_bound_coef,
  estimate_upper = 1,
  upper = 1
)

results[["LATE"]] <- tibble(
  estimand = "LATE",
  lower = lpr_lb_LATE,
  estimate_lower = lpr2_coef,
  estimate_upper = lpr2_coef,
  upper = min(1, lpr2_coef + cv_cns2 * lpr2_se)
)

results[["LATE YZ"]] <- tibble(
  estimand = "LATE: Y,Z observed",
  lower = max(0, itt_y["estimate"] - cv_cns1 * itt_y["se"]),
  estimate_lower = itt_y["estimate"],
  estimate_upper = 1,
  upper = 1
)

# APR conditional on voting without persuasive treatment.
# Stata uses voted, which is expected to exist in publicdata.dta.
if ("voted" %in% names(gkb)) {

  gkb <- gkb %>%
    mutate(yt00 = (1 - voted) * (1 - treat))

  multi_point <- function(data) {
    fit_num <- lm(outcome ~ instr, data = data)
    fit_den <- lm(yt00 ~ instr, data = data)
    b_num <- coef(fit_num)
    b_den <- coef(fit_den)
    b_num[["instr"]] /
      (1 - b_num[["(Intercept)"]] - b_den[["(Intercept)"]])
  }

  multi_lower_bound_coef <- multi_point(gkb)

  set.seed(12345)
  multi_boot <- replicate(
    500,
    {
      idx <- sample.int(nrow(gkb), replace = TRUE)
      multi_point(gkb[idx, ])
    }
  )

  multi_lower_bound_se <- sd(multi_boot, na.rm = TRUE)

  results[["APR conditional"]] <- tibble(
    estimand = "APR conditional on voting without treatment",
    lower = max(0, multi_lower_bound_coef - cv_cns1_adj * multi_lower_bound_se),
    estimate_lower = multi_lower_bound_coef,
    estimate_upper = 1,
    upper = 1
  )
} else {
  # Uncertain revision: the trace references variable voted. If the R copy of
  # publicdata.dta does not contain voted, this row is omitted rather than
  # silently substituting another variable.
  warning("Variable `voted` not found; omitted conditional APR row.")
}

table4_data <- bind_rows(results) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

table4 <- tt(
  table4_data,
  caption = "Table 4: Persuasion Rate Bounds and Related Estimands"
)

save_tt(
  table4,
  file = file.path(results_dir, "table4.tex")
)

write.csv(
  table4_data,
  file.path(results_dir, "table4.csv"),
  row.names = FALSE
)

print(table4)
