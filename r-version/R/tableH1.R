# ************
# * SCRIPT:   tableH1.do
# * PURPOSE:  Creates Table H1
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
# * SCRIPT:   tableH1.R
# * PURPOSE:  Creates Table H1
# ************

# Original Stata source: trace_tableH1.txt.
# Translates survey-weighted probit workflow into weighted glm(..., probit).
# The final SE uses a weighted influence-function mean approximation.
# This is the main uncertain revision because Stata's svy: mean linearization
# and R's base glm do not share a single native drop-in equivalent.

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

weighted_quantile <- function(x, w, prob = 0.5) {
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= prob)[1]]
}

weighted_mean <- function(x, w) {
  weighted.mean(x, w, na.rm = TRUE)
}

weighted_mean_se <- function(x, w) {
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]
  mu <- weighted.mean(x, w)
  sqrt(sum(w^2 * (x - mu)^2) / sum(w)^2)
}

ntv <- ntv %>%
  mutate(
    y_vote_Unity = as.integer(vote_Unity == 0),
    y_vote_OVR   = as.integer(vote_OVR == 1)
  )

IV_p50 <- weighted_quantile(
  ntv$tvmaxtveloss5050powerA,
  ntv$kishweig,
  prob = 0.5
)

ntv <- ntv %>%
  mutate(
    IV_NTV = as.integer(tvmaxtveloss5050powerA > IV_p50)
  )

z_formula <- as.formula(
  paste(
    "IV_NTV ~",
    paste(c(sociodem, basic), collapse = " + ")
  )
)

z_fit <- glm(
  z_formula,
  data = ntv,
  weights = kishweig,
  family = binomial(link = "probit")
)

ntv$phat_Z1 <- predict(
  z_fit,
  newdata = ntv,
  type = "response"
)

t_formula <- as.formula(
  paste(
    "Watches_NTV_1999 ~",
    paste(c(sociodem, basic), collapse = " + ")
  )
)

t_fit_z1 <- glm(
  t_formula,
  data = ntv %>% filter(IV_NTV == 1),
  weights = kishweig,
  family = binomial(link = "probit")
)

ntv$phat_T_Z1 <- predict(
  t_fit_z1,
  newdata = ntv,
  type = "response"
)

t_fit_z0 <- glm(
  t_formula,
  data = ntv %>% filter(IV_NTV == 0),
  weights = kishweig,
  family = binomial(link = "probit")
)

ntv$phat_T_Z0 <- predict(
  t_fit_z0,
  newdata = ntv,
  type = "response"
)

tableH1_results <- list()

for (party in party_list) {

  y_vote <- paste0("y_vote_", party)
  newvote <- paste0("newvote_", party)

  ntv[[newvote]] <- as.integer(ntv[[y_vote]] == 1)

  y_formula <- as.formula(
    paste(
      newvote, "~",
      paste(c(sociodem, basic), collapse = " + ")
    )
  )

  y_fit_z1 <- glm(
    y_formula,
    data = ntv %>% filter(IV_NTV == 1),
    weights = kishweig,
    family = binomial(link = "probit")
  )

  ntv[[paste0("phat_Y_Z1_", party)]] <- predict(
    y_fit_z1,
    newdata = ntv,
    type = "response"
  )

  y_fit_z0 <- glm(
    y_formula,
    data = ntv %>% filter(IV_NTV == 0),
    weights = kishweig,
    family = binomial(link = "probit")
  )

  ntv[[paste0("phat_Y_Z0_", party)]] <- predict(
    y_fit_z0,
    newdata = ntv,
    type = "response"
  )

  phat_y_z1 <- ntv[[paste0("phat_Y_Z1_", party)]]
  phat_y_z0 <- ntv[[paste0("phat_Y_Z0_", party)]]

  w1 <- ntv$IV_NTV / ntv$phat_Z1
  w0 <- (1 - ntv$IV_NTV) / (1 - ntv$phat_Z1)

  t1 <- w1 * (ntv[[newvote]] - phat_y_z1)
  t0 <- w0 * (ntv[[newvote]] - phat_y_z0)

  beta1 <- weighted_mean(phat_y_z1, ntv$kishweig)
  beta0 <- weighted_mean(phat_y_z0, ntv$kishweig)

  ate_lb <- (beta1 - beta0) / (1 - beta0)

  infl_t1 <- t1 + (phat_y_z1 - beta1)
  infl_t0 <- t0 + (phat_y_z0 - beta0)

  infl <- (1 / (1 - beta0)) *
    (
      infl_t1 -
        ((1 - beta1) / (1 - beta0)) * infl_t0
    )

  se <- weighted_mean_se(infl, ntv$kishweig)

  tableH1_results[[party]] <- tibble(
    Party = party,
    Estimate = ate_lb,
    `Std. Error` = se,
    `80% CI Lower` = ate_lb - 1.645 * se
  )
}

tableH1_data <- bind_rows(tableH1_results) %>%
  mutate(across(where(is.numeric), ~ round(100 * .x, 2)))

tableH1 <- tt(
  tableH1_data,
  caption = "Average Persuasion Rate with Binary Instrument"
)

save_tt(
  tableH1,
  file = file.path(results_dir, "tableH1.tex")
)

write.csv(
  tableH1_data,
  file.path(results_dir, "tableH1.csv"),
  row.names = FALSE
)

print(tableH1)
