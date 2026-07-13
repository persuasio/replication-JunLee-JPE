# ************
# * SCRIPT:   tableD1.do
# * PURPOSE:  Creates Table D1
# *
# * ACKNOWLEDGMENT
# *       The orginal dataset "FoxNewsFinalDataQJE" is from DellaVigna and Kaplan (QJE, 2007).
# *       The dataset is available at http://eml.berkeley.edu/~sdellavi/index.html.
# *   We thank the authors of the paper to make their data available online.
# *
# * NOTES
# *       An alternative measure of the dependent variable is used.
# *   That is, Y = (Republican voting share as a share of the voting-age population)
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
# * SCRIPT:   tableD1.R
# * PURPOSE:  Creates Table D1
# ************

library(fixest)

# Original Stata source: trace_tableD1.txt.
# Uses DellaVigna and Kaplan data and fixed-effect regressions to construct
# persuasion-rate bounds. Stata areg is translated to fixest::feols.

fox <- read_dta(
  file.path(
    persuasion_dir,
    "data",
    "DellaVignaKaplan2007",
    "FoxNewsFinalDataQJE.dta"
  )
)

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

# Stata shorthand: poptot2000d2-poptot2000d10 noch2000d2-noch2000d10.
contrcbl2000 <- c(
  paste0("poptot2000d", 2:10),
  paste0("noch2000d", 2:10)
)

fox <- fox %>%
  filter(sample12000 == 1) %>%
  mutate(
    new_y2000 = reppresvotes2000 / pop18p2000,
    new_y1996 = reppresvotes1996 / pop18p1996
  )

rhs <- c(
  "new_y1996",
  "foxnews2000",
  contrcens2000,
  contrcens00m90,
  contrcbl2000
)

results <- list()

for (fixed_effect in c("diststate", "countystate")) {

  # Stata: areg new_y2000 ... , a(fixed_effect) robust cluster(account2000)
  y_formula <- as.formula(
    paste(
      "new_y2000 ~",
      paste(rhs, collapse = " + "),
      "|",
      fixed_effect
    )
  )

  y_fit <- feols(
    y_formula,
    data = fox,
    weights = ~ pop18p1996,
    cluster = ~ account2000
  )

  fox$yhat <- predict(y_fit, newdata = fox)

  coeff_b <- coef(y_fit)[["foxnews2000"]]

  fox <- fox %>%
    mutate(
      yhat1 = clip01(yhat + coeff_b - coeff_b * foxnews2000),
      yhat0 = clip01(yhat - coeff_b * foxnews2000),
      thetahat_num = yhat1 - yhat0,
      thetahat_den = pmax(1 - yhat0, 1e-8)
    )

  avg_num <- weighted.mean(
    fox$thetahat_num[fox$auddiaryScar > 0],
    fox$auddiaryScar[fox$auddiaryScar > 0],
    na.rm = TRUE
  )

  avg_den <- weighted.mean(
    fox$thetahat_den[fox$auddiaryScar > 0],
    fox$auddiaryScar[fox$auddiaryScar > 0],
    na.rm = TRUE
  )

  avg_lb <- avg_num / avg_den

  # Stata: areg foxanyScar foxnews2000 ...
  t_formula <- as.formula(
    paste(
      "foxanyScar ~",
      paste(c("foxnews2000", contrcens2000, contrcens00m90, contrcbl2000), collapse = " + "),
      "|",
      fixed_effect
    )
  )

  t_fit <- feols(
    t_formula,
    data = fox,
    weights = ~ auddiaryScar,
    cluster = ~ account2000
  )

  fox$ehat <- predict(t_fit, newdata = fox)

  coeff_e <- coef(t_fit)[["foxnews2000"]]

  fox <- fox %>%
    mutate(
      ehat1 = clip01(ehat + coeff_e - coeff_e * foxnews2000),
      ehat0 = clip01(ehat - coeff_e * foxnews2000),
      ub_num1 = yhat1 + 1 - ehat1,
      ub_num2 = yhat0 - ehat0,
      ub_num = pmin(1, ub_num1) - pmax(0, ub_num2),
      ub_den = pmax(1 - pmax(0, ub_num2), 1e-8),
      late_weight = ehat1 - ehat0,
      theta_local_den = pmin(thetahat_den, late_weight)
    )

  ub_avg_num <- weighted.mean(
    fox$ub_num[fox$auddiaryScar > 0],
    fox$auddiaryScar[fox$auddiaryScar > 0],
    na.rm = TRUE
  )

  ub_avg_den <- weighted.mean(
    fox$ub_den[fox$auddiaryScar > 0],
    fox$auddiaryScar[fox$auddiaryScar > 0],
    na.rm = TRUE
  )

  local_den <- weighted.mean(
    fox$theta_local_den[fox$auddiaryScar > 0],
    fox$auddiaryScar[fox$auddiaryScar > 0],
    na.rm = TRUE
  )

  results[[fixed_effect]] <- tibble(
    fixed_effect = fixed_effect,
    `APR (LB)` = avg_lb,
    `APR (UB)` = ub_avg_num / ub_avg_den,
    `LPR (LB)` = avg_num / local_den,
    `LPR (UB)` = 1
  )

  fox <- fox %>%
    select(
      -yhat, -yhat1, -yhat0, -thetahat_num, -thetahat_den,
      -ehat, -ehat1, -ehat0, -starts_with("ub_"),
      -late_weight, -theta_local_den
    )
}

tableD1_data <- bind_rows(results) %>%
  mutate(across(where(is.numeric), ~ round(100 * .x, 2)))

tableD1 <- tt(
  tableD1_data,
  caption = "Presidential Effects: Alternative Dependent Variable"
)

save_tt(
  tableD1,
  file = file.path(results_dir, "tableD1.tex")
)

write.csv(
  tableD1_data,
  file.path(results_dir, "tableD1.csv"),
  row.names = FALSE
)

print(tableD1)
