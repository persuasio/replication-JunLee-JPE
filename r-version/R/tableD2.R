# ************
# * SCRIPT:   tableD2.do
# * PURPOSE:  Creates Table D2
# *
# * ACKNOWLEDGMENT
# *       The orginal dataset "NTV_Individual_Data.dta" is from Enikolopov, Petrova, and Zhuravskaya (AER, 2011).
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
# * SCRIPT:   tableD2.R
# * PURPOSE:  Creates Table D2
# ************

library(writexl)

# Original Stata source: trace_tableD2.txt.
# This is the same marginal persuasion-rate workflow as figureD1input.R, with
# the final matrix exported as Table D2.

ntv <- read_dta(
  file.path(
    persuasion_dir,
    "data",
    "EnikolopovPetrovaZhuravskaya2011",
    "NTV_Individual_Data.dta"
  )
)

basic <- c("logpop98", "wage98_ln")
sociodem <- c("male", "age", "educ1", "married", "consump")
party_list <- c("Unity", "OVR")

ntv <- ntv %>%
  mutate(
    y_vote_Unity = as.integer(vote_Unity == 0),
    y_vote_OVR   = as.integer(vote_OVR == 1)
  )

exposure_fit <- lm(
  Watches_NTV_1999 ~
    male + age + educ1 + married + consump +
    logpop98 + wage98_ln +
    tvmaxtveloss5050powerA +
    I(tvmaxtveloss5050powerA^2) +
    I(tvmaxtveloss5050powerA^3) +
    tvmaxtveloss5050powerA:male +
    tvmaxtveloss5050powerA:age +
    tvmaxtveloss5050powerA:educ1 +
    tvmaxtveloss5050powerA:married +
    tvmaxtveloss5050powerA:consump,
  data = ntv,
  weights = kishweig
)

ntv$phat <- clip01(predict(exposure_fit, newdata = ntv))

mte_grid <- seq(0.40, 0.60, by = 0.01)

marginal_effect_at_phat <- function(fit, data) {

  b <- coef(fit)

  vapply(mte_grid, function(v) {

    d <- rep(0, nrow(data))

    d <- d + ifelse("phat" %in% names(b), b[["phat"]], 0)
    d <- d + ifelse("I(phat^2)" %in% names(b), 2 * b[["I(phat^2)"]] * v, 0)
    d <- d + ifelse("I(phat^3)" %in% names(b), 3 * b[["I(phat^3)"]] * v^2, 0)

    for (x in c("male", "age", "educ1", "married", "consump")) {
      term1 <- paste0("phat:", x)
      term2 <- paste0(x, ":phat")
      bx <- if (term1 %in% names(b)) {
        b[[term1]]
      } else if (term2 %in% names(b)) {
        b[[term2]]
      } else {
        0
      }
      d <- d + bx * data[[x]]
    }

    weighted.mean(d, data$kishweig, na.rm = TRUE)
  }, numeric(1))
}

results_long <- vector("list", length(party_list))

for (party in party_list) {

  y_var <- paste0("y_vote_", party)
  notwatch_var <- paste0("notwatch_vote_", party)

  ntv[[notwatch_var]] <- as.integer(
    ntv[[y_var]] == 1 &
      ntv$Watches_NTV_1999 == 0
  )

  y_fit <- lm(
    as.formula(
      paste0(
        y_var,
        " ~ male + age + educ1 + married + consump + logpop98 + wage98_ln + ",
        "phat + I(phat^2) + I(phat^3) + ",
        "phat:male + phat:age + phat:educ1 + phat:married + phat:consump"
      )
    ),
    data = ntv,
    weights = kishweig
  )

  yw_fit <- lm(
    as.formula(
      paste0(
        notwatch_var,
        " ~ male + age + educ1 + married + consump + logpop98 + wage98_ln + ",
        "phat + I(phat^2) + I(phat^3) + ",
        "phat:male + phat:age + phat:educ1 + phat:married + phat:consump"
      )
    ),
    data = ntv,
    weights = kishweig
  )

  num <- marginal_effect_at_phat(y_fit, ntv)
  den <- 1 + marginal_effect_at_phat(yw_fit, ntv)

  results_long[[party]] <- tibble(
    row = sprintf("v = %.2f", mte_grid),
    party = party,
    estimate = num / den
  )
}

tableD2 <- bind_rows(results_long) %>%
  bind_rows(
    bind_rows(results_long) %>%
      group_by(party) %>%
      summarise(
        row = "Avg between 0.4 and 0.6",
        estimate = mean(estimate, na.rm = TRUE),
        .groups = "drop"
      )
  ) %>%
  tidyr::pivot_wider(
    names_from = party,
    values_from = estimate
  ) %>%
  arrange(
    ifelse(row == "Avg between 0.4 and 0.6", 0, 1),
    row
  )

write.csv(
  tableD2,
  file.path(results_dir, "tableD2.csv"),
  row.names = FALSE
)

writexl::write_xlsx(
  tableD2,
  file.path(results_dir, "tableD2.xlsx")
)

tableD2_tbl <- tt(
  tableD2,
  caption = "Estimates of Marginal and Average Persuasion Rates"
)

save_tt(
  tableD2_tbl,
  file = file.path(results_dir, "tableD2.tex")
)

print(tableD2_tbl)
