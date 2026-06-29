# ************
# * SCRIPT:   tableE2.do
# * PURPOSE:  Creates Table E2
# *
# * ACKNOWLEDGMENT
# *       The orginal dataset "CharityOutputQJE" is from DellaVigna, List, and Malmendier (QJE, 2012).
# *       The data collected in DLM are available at http://eml.berkeley.edu/~sdellavi/index.html
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
# * SCRIPT:   tableE2.R
# * PURPOSE:  Creates Table E2
# ************

# Original Stata source: trace_tableE2.txt.
# This script uses fitted values from linear probability models with fixed-effect
# dummies. Stata xi: reg is translated to lm(... + factor(...)).
# Cluster-robust SEs in Stata are not used in the final table calculations here,
# which are based on fitted-value averages.

charity <- read_dta(
  file.path(
    persuasion_dir,
    "data",
    "DellaVignaListMalmendier2012",
    "CharityOutputQJE.dta"
  )
)

charity <- charity %>%
  mutate(treatment = ifelse(treatment == "2Ww", "W", treatment)) %>%
  filter(nosol != 1) %>%
  select(-nosol) %>%
  filter(toeliminate != 1) %>%
  select(-toeliminate) %>%
  filter(!(sol == "Angelena" & date %in% as.Date(c("2008-07-27", "2008-07-13")))) %>%
  filter(!(sol == "Shedora"  & date %in% as.Date(c("2008-07-27", "2008-08-10")))) %>%
  filter(!(sol == "Tehmur"   & date == as.Date("2008-06-01"))) %>%
  filter(!(sol == "Phillip"  & date == as.Date("2008-08-09"))) %>%
  filter(!(sol == "Robert"   & date == as.Date("2008-07-13") & hour %in% c(11, 13))) %>%
  mutate(
    month = as.integer(format(date, "%m")),
    dwave = case_when(
      (month %in% c(7, 8) & charity == "LaRabida") & year == 2008 ~ 1,
      ((month %in% c(4, 5, 6)) | (charity == "Ecu" & month == 7)) & year == 2008 ~ 0,
      (month %in% c(9, 10)) & year == 2008 ~ 2,
      TRUE ~ NA_real_
    ),
    sodate = paste(solicitor, date, sep = "_"),
    treatmentby = paste(treatment, charity, sep = "_"),
    grsol = as.integer(factor(solicitor)),
    grdatloc = as.integer(factor(paste(date, location))),
    grdatlocsol = as.integer(factor(paste(date, location, solicitor))),
    grhour = as.integer(factor(hour)),
    grarea = as.integer(factor(area_rank)),
    amt_donate = ifelse(is.na(amt_donate), 0, amt_donate),
    dW = as.integer(treatment == "W"),
    dOo = as.integer(treatment == "Oo"),
    dWEcu = as.integer(treatment == "W" & charity == "Ecu"),
    dWLar = as.integer(treatment == "W" & charity == "LaRabida"),
    dOoEcu = as.integer(treatment == "Oo" & charity == "Ecu"),
    dOoLar = as.integer(treatment == "Oo" & charity == "LaRabida"),
    dEcu = as.integer(charity == "Ecu")
  )

for (x in c("0d5m", "0d10m", "5d5m", "5d10m", "10d10m", "10d5m")) {
  for (y in c("Nw", "W", "Oo")) {
    charity$treatment[charity$treatment == paste0(y, "-", x)] <- paste0(y, x)
  }
}

for (x in c("Nw0d10m", "W0d10m", "W0d5m", "W10d10m")) {
  charity[[paste0("d", x, "08")]] <- as.integer(
    charity$treatment == x &
      charity$year == 2008
  )
}

for (x in c("Nw0d5m", "Nw5d5m", "W0d10m", "W0d5m", "W10d5m", "W5d5m", "Oo0d5m", "Oo5d5m")) {
  charity[[paste0("d", x, "09")]] <- as.integer(
    charity$treatment == x &
      charity$year == 2009
  )
}

# Stata drops omitted categories dNw0d10m08 and dNw0d5m09.
charity <- charity %>%
  select(-any_of(c("dNw0d10m08", "dNw0d5m09")))

charity_results <- list()

for (ch in c("LaRabida", "Ecu")) {

  tmp <- charity %>%
    filter(charity == ch)

  a_fit <- lm(
    answer ~ dW + dOo + factor(grsol) + factor(grdatloc) + factor(grhour) + factor(grarea),
    data = tmp
  )

  tmp$ahat <- predict(a_fit, newdata = tmp)
  coeff_a_W <- coef(a_fit)[["dW"]]
  coeff_a_Oo <- coef(a_fit)[["dOo"]]

  tmp <- tmp %>%
    mutate(
      ahat_Nw = clip01(ahat - coeff_a_W * dW - coeff_a_Oo * dOo),
      ahat_W  = clip01(ahat_Nw + coeff_a_W),
      ahat_Oo = clip01(ahat_Nw + coeff_a_Oo)
    )

  y_fit <- lm(
    saidyes ~ dW + dOo + factor(grsol) + factor(grdatloc) + factor(grhour) + factor(grarea),
    data = tmp
  )

  tmp$yhat <- predict(y_fit, newdata = tmp)
  coeff_y_W <- coef(y_fit)[["dW"]]
  coeff_y_Oo <- coef(y_fit)[["dOo"]]

  tmp <- tmp %>%
    mutate(
      yhat_Nw = clip01(yhat - coeff_y_W * dW - coeff_y_Oo * dOo),
      yhat_W  = clip01(yhat_Nw + coeff_y_W),
      yhat_Oo = clip01(yhat_Nw + coeff_y_Oo),
      lb_Nw = yhat_Nw,
      lb_W  = yhat_W,
      lb_Oo = yhat_Oo,
      ub_Nw = yhat_Nw + 1 - ahat_Nw,
      ub_W  = yhat_W  + 1 - ahat_W,
      ub_Oo = yhat_Oo + 1 - ahat_Oo
    )

  for (arm in c("Nw", "W", "Oo")) {

    yhat_mean <- mean(tmp[[paste0("yhat_", arm)]], na.rm = TRUE)
    ahat_mean <- mean(tmp[[paste0("ahat_", arm)]], na.rm = TRUE)
    lb_mean   <- mean(tmp[[paste0("lb_", arm)]], na.rm = TRUE)
    ub_mean   <- mean(tmp[[paste0("ub_", arm)]], na.rm = TRUE)

    charity_results[[paste(ch, arm)]] <- tibble(
      Charity = ch,
      Treatment = arm,
      `P(Y=1|Z=1)` = yhat_mean,
      `e(1)` = ahat_mean,
      `APR (LB)` = lb_mean,
      `APR (UB)` = ub_mean,
      LPR = yhat_mean / ahat_mean
    )
  }
}

tableE2_data <- bind_rows(charity_results) %>%
  mutate(across(where(is.numeric), ~ round(100 * .x, 1)))

tableE2 <- tt(
  tableE2_data,
  caption = "Persuasive Effect by Treatment in DellaVigna, List, and Malmendier (2012)"
)

save_tt(
  tableE2,
  file = file.path(results_dir, "tableE2.tex")
)

write.csv(
  tableE2_data,
  file.path(results_dir, "tableE2.csv"),
  row.names = FALSE
)

print(tableE2)
