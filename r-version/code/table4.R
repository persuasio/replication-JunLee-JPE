# ************
# * SCRIPT:   table4.R
# * PURPOSE:  Creates Table 4
# *
# * ACKNOWLEDGMENT
# *   The original dataset "publicdata.dta" is from
# *   Gerber, Karlan, and Bergan (2009, AEJ Applied).
# ************

library(haven)
library(dplyr)
library(tibble)
library(tinytable)
library(here)
library(sandwich)
library(lmtest)
library(AER)
library(numDeriv)

output_dir <- here::here("output")
data_path <- here::here(
  "data",
  "GerberKarlanBergan2009",
  "publicdata.dta"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ",
    data_path,
    "\nRun this script from within the replication project, ",
    "with data/GerberKarlanBergan2009/publicdata.dta present."
  )
}

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

# The Stata code manually constructs ITT, APR, LPR, and LATE bounds.
# This R version follows those formulas directly. For suest/nlcom sections,
# we use numerical delta-method SEs with HC1 covariance matrices.

gkb <- read_dta(data_path) %>%
  filter(
    .data$times != 1,
    !is.na(.data$voteddem_all),
    !is.na(.data$readsome)
  ) %>%
  mutate(
    case_id = row_number(),
    outcome = .data$voteddem_all,
    treat = .data$readsome,
    instr = .data$post,
    a_u = .data$outcome * .data$treat + (1 - .data$treat),
    b_l = .data$outcome * (1 - .data$treat),
    a_u2 = .data$outcome + (1 - .data$treat),
    b_l2 = .data$outcome - .data$treat
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

# ITT estimate on outcome
itt_y_fit <- lm(outcome ~ instr, data = gkb)
itt_y <- robust_coef(itt_y_fit, "instr")

results[["ITT outcome"]] <- tibble(
  estimand = "ITT: outcome",
  lower = max(0, itt_y["estimate"] - cv_cns2 * itt_y["se"]),
  estimate_lower = itt_y["estimate"],
  estimate_upper = itt_y["estimate"],
  upper = min(1, itt_y["estimate"] + cv_cns2 * itt_y["se"])
)

# ITT estimate on treatment
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
gkb <- gkb %>%
  mutate(yt00 = (1 - .data$voted) * (1 - .data$treat))

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

table4_results <- bind_rows(results)


stopifnot(nrow(table4_results) == 11L)

# Preserve all 11 numeric results for verification and downstream use.
table4_raw <- table4_results %>%
  transmute(
    Estimand = estimand,
    `Estimate (LB)` = estimate_lower,
    `Estimate (UB)` = estimate_upper,
    `CI (LB)` = lower,
    `CI (UB)` = upper
  )

# Reformat the table like the paper
fmt_num <- function(x) {
  x <- as.numeric(x)

  if (is.na(x)) {
    return("")
  }

  if (isTRUE(all.equal(x, 0))) {
    return("0")
  }

  if (isTRUE(all.equal(x, 1))) {
    return("1")
  }

  out <- sprintf("%.4f", x)
  out <- sub("^0\\.", ".", out)
  out <- sub("^-0\\.", "-.", out)
  out
}

fmt_interval <- function(lb, ub) {
  paste0("[", fmt_num(lb), ", ", fmt_num(ub), "]")
}

fmt_ci <- function(lb, ub) {
  paste0("{", fmt_num(lb), "–", fmt_num(ub), "}")
}

fmt_point_cell <- function(est, lb, ub) {
  paste0(
    fmt_num(est),
    "<br>",
    fmt_ci(lb, ub)
  )
}

fmt_interval_cell <- function(est_lb, est_ub, ci_lb, ci_ub) {
  paste0(
    fmt_interval(est_lb, est_ub),
    "<br>",
    fmt_ci(ci_lb, ci_ub)
  )
}

get_result <- function(name) {
  row <- table4_results %>%
    filter(.data$estimand == name)

  stopifnot(nrow(row) == 1L)
  row
}

itt <- get_result("ITT: outcome")

avg_1 <- get_result("APR: Y,T,Z observed")
avg_2 <- get_result("APR: Y,T and Y,Z observed")
avg_3 <- get_result("APR: Y,Z observed")

late_1 <- get_result("LATE")
late_3 <- get_result("LATE: Y,Z observed")

local_1 <- get_result("LPR: Y,T,Z observed")
local_2 <- get_result("LPR: Y,T and Y,Z observed")
local_3 <- get_result("LPR: Y,Z observed")

# Publication-facing Table 4:
# estimates or identified sets first, confidence intervals underneath.
table4_display <- tibble(
  Parameter = c(
    "ITT",
    "\u03b8_avg",
    "LATE",
    "\u03b8_local"
  ),

  `(Y_i, T_i, Z_i)` = c(
    "",
    fmt_interval_cell(
      avg_1$estimate_lower,
      avg_1$estimate_upper,
      avg_1$lower,
      avg_1$upper
    ),
    fmt_point_cell(
      late_1$estimate_lower,
      late_1$lower,
      late_1$upper
    ),
    fmt_point_cell(
      local_1$estimate_lower,
      local_1$lower,
      local_1$upper
    )
  ),

  `(Y_i, Z_i) and (T_i, Z_i)` = c(
    fmt_point_cell(
      itt$estimate_lower,
      itt$lower,
      itt$upper
    ),
    fmt_interval_cell(
      avg_2$estimate_lower,
      avg_2$estimate_upper,
      avg_2$lower,
      avg_2$upper
    ),
    "",
    fmt_interval_cell(
      local_2$estimate_lower,
      local_2$estimate_upper,
      local_2$lower,
      local_2$upper
    )
  ),

  `(Y_i, Z_i) Only` = c(
    "",
    fmt_interval_cell(
      avg_3$estimate_lower,
      avg_3$estimate_upper,
      avg_3$lower,
      avg_3$upper
    ),
    fmt_interval_cell(
      late_3$estimate_lower,
      late_3$estimate_upper,
      late_3$lower,
      late_3$upper
    ),
    fmt_interval_cell(
      local_3$estimate_lower,
      local_3$estimate_upper,
      local_3$lower,
      local_3$upper
    )
  )
)

write.csv(
  table4_raw,
  file = file.path(output_dir, "table4_raw.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  table4_display,
  file = file.path(output_dir, "table4_display.csv"),
  row.names = FALSE,
  na = ""
)

table4_tex <- tt(
  table4_display,
  caption = "Estimates of Key Parameters"
)

print(table4_tex)

save_tt(
  table4_tex,
  output = file.path(output_dir, "table4.tex"),
  overwrite = TRUE
)
