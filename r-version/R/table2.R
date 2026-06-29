# ************
# * SCRIPT:   table2.R
# * PURPOSE:  Creates Table 2
# ************

library(dplyr)
library(tibble)
library(broom)
library(lmtest)
library(sandwich)
library(persuasio)
library(tinytable)

# This script follows the Stata table2.do workflow:
#   1. estimate the complier share from active_user on treatment assignment;
#   2. compute the CY19 persuasion rate from the reduced-form coefficient;
#   3. compute APR and LPR using the R persuasio package;
#   4. reorder rows and multiply by 100 before exporting.
#
# It assumes the cleaned CY19 dataset has already been created by CY19_data.R

if (!exists("ChenYang2019.rds")) {
  stop("Object was not found. Run/source CY19_data.R first.")
}
cy <- "ChenYang2019.rds"

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# Number of bootstrap replications. The Stata trace used 50 because $nbt was
# empty in that run. Increase for final output if desired.
nbt <- 50

# -----------------------------------------------------------------------------
# Helper: extract estimates from persuasio() output
# -----------------------------------------------------------------------------

extract_persuasio_estimates <- function(out, estimand) {

  # persuasio objects may print nicely while storing estimates in slightly
  # different list components across versions. This extractor searches common
  # locations without relying on the Stata trace.

  candidates <- list(
    out$estimates,
    out$estimate,
    out$results,
    out$result,
    out$table,
    out
  )

  for (x in candidates) {

    if (is.null(x)) {
      next
    }

    df <- tryCatch(
      as.data.frame(x),
      error = function(e) NULL
    )

    if (is.null(df) || nrow(df) == 0) {
      next
    }

    names_lower <- tolower(names(df))

    if (estimand == "apr") {

      lb_col <- which(names_lower %in% c(
        "lower bound", "lower_bound", "apr_lb", "lb", "lower"
      ))

      ub_col <- which(names_lower %in% c(
        "upper bound", "upper_bound", "apr_ub", "ub", "upper"
      ))

      cil_col <- which(names_lower %in% c(
        "ci lower", "ci_lower", "ci_lb", "conf.low", "conf_low"
      ))

      ciu_col <- which(names_lower %in% c(
        "ci upper", "ci_upper", "ci_ub", "conf.high", "conf_high"
      ))

      if (length(lb_col) > 0 && length(ub_col) > 0) {
        return(
          tibble(
            apr_lb = as.numeric(df[[lb_col[1]]][1]),
            apr_ub = as.numeric(df[[ub_col[1]]][1]),
            apr_ci_lb = if (length(cil_col) > 0) as.numeric(df[[cil_col[1]]][1]) else NA_real_,
            apr_ci_ub = if (length(ciu_col) > 0) as.numeric(df[[ciu_col[1]]][1]) else NA_real_
          )
        )
      }
    }

    if (estimand == "lpr") {

      lpr_col <- which(names_lower %in% c(
        "lpr", "estimate", "local persuasion rate", "local_persuasion_rate"
      ))

      cil_col <- which(names_lower %in% c(
        "ci lower", "ci_lower", "ci_lb", "conf.low", "conf_low"
      ))

      ciu_col <- which(names_lower %in% c(
        "ci upper", "ci_upper", "ci_ub", "conf.high", "conf_high"
      ))

      if (length(lpr_col) > 0) {
        return(
          tibble(
            lpr = as.numeric(df[[lpr_col[1]]][1]),
            lpr_ci_lb = if (length(cil_col) > 0) as.numeric(df[[cil_col[1]]][1]) else NA_real_,
            lpr_ci_ub = if (length(ciu_col) > 0) as.numeric(df[[ciu_col[1]]][1]) else NA_real_
          )
        )
      }
    }
  }

  stop(
    paste0(
      "Could not extract ", estimand,
      " estimates from persuasio() output. Inspect str(out) and update extract_persuasio_estimates()."
    )
  )
}

# -----------------------------------------------------------------------------
# Helper: one Table 2 row
# -----------------------------------------------------------------------------

make_table2_row <- function(data, y, y0 = NULL) {

  # Stata: qui reg active_user treatment_vpnonly treatment_nlonly treatment_vpnnl, r
  active_fit <- lm(
    active_user ~ treatment_vpnonly + treatment_nlonly + treatment_vpnnl,
    data = data
  )

  share_compliers <- coef(active_fit)[["treatment_vpnnl"]]

  # Stata panel variables use Y_w3_p as outcome and Y_w1_p among treatment_vpnnl
  # as the pre-treatment baseline share. Non-panel variables use the intercept.
  rf_fit <- lm(
    reformulate(
      c("treatment_vpnonly", "treatment_nlonly", "treatment_vpnnl"),
      response = y
    ),
    data = data
  )

  b_vpnnl <- coef(rf_fit)[["treatment_vpnnl"]]

  if (!is.null(y0)) {
    baseline_share <- mean(
      data[[y0]][data$treatment_vpnnl == 1],
      na.rm = TRUE
    )
  } else {
    baseline_share <- coef(rf_fit)[["(Intercept)"]]
  }

  persuasion_DK <- b_vpnnl / (share_compliers * (1 - baseline_share))

  # Stata: persuasio apr y active_user treatment_vpnnl if treatment_main == 1 | treatment_main == 3
  # R: use the eponymous R package API from README.md.
  pers_data <- data %>%
    filter(
      treatment_main == 1 |
        treatment_main == 3
    )

  apr_out <- persuasio(
    est = "apr",
    y = y,
    t = "active_user",
    z = "treatment_vpnnl",
    data = pers_data,
    method = "bootstrap",
    nboot = nbt
  )

  lpr_out <- persuasio(
    est = "lpr",
    y = y,
    t = "active_user",
    z = "treatment_vpnnl",
    data = pers_data,
    method = "bootstrap",
    nboot = nbt
  )

  apr <- extract_persuasio_estimates(apr_out, "apr")
  lpr <- extract_persuasio_estimates(lpr_out, "lpr")

  tibble(
    variable = y,
    DK = persuasion_DK,
    `APR (LB)` = apr$apr_lb,
    `APR (UB)` = apr$apr_ub,
    LPR = lpr$lpr,
    LPR_CI_LB = lpr$lpr_ci_lb,
    LPR_CI_UB = lpr$lpr_ci_ub,
    APR_CI_LB = apr$apr_ci_lb,
    APR_CI_UB = apr$apr_ci_ub
  )
}

# -----------------------------------------------------------------------------
# Stata Table 2 variable order
# -----------------------------------------------------------------------------

panel_vars <- c(
  "info_foreign_website",
  "info_freq_website_for",
  "az_belief_media_value",
  "az_belief_media_trust",
  "bias_domestic",
  "bias_foreign",
  "az_belief_media_justif",
  "bias_dom_govt_policy_t1",
  "bias_for_govt_policy_t1"
)

nonpanel_vars <- c(
  "vpn_purchase_wmt_record",
  "vpn_purchase_yes"
)

panel_results <- lapply(panel_vars, function(Y) {
  make_table2_row(
    data = cy,
    y = paste0(Y, "_w3_p"),
    y0 = paste0(Y, "_w1_p")
  )
})

nonpanel_results <- lapply(nonpanel_vars, function(Y) {
  make_table2_row(
    data = cy,
    y = paste0(Y, "_p"),
    y0 = NULL
  )
})

results <- bind_rows(
  panel_results,
  nonpanel_results
)

# first four panel variables, then the two purchase variables, then remaining panel variables.
results <- bind_rows(
  results[1:4, ],
  results[10:11, ],
  results[5:9, ]
)

# Stata: matrix results = 100 * results
results_pct <- results %>%
  mutate(
    across(
      c(DK, `APR (LB)`, `APR (UB)`, LPR, LPR_CI_LB, LPR_CI_UB, APR_CI_LB, APR_CI_UB),
      ~ 100 * .x
    )
  )

# Main printed table follows the Stata frmttable columns:
# CY19, APR (LB), APR (UB), LPR, plus LPR CI columns retained explicitly.
table2 <- results_pct %>%
  transmute(
    Outcome = variable,
    CY19 = round(DK, 1),
    `APR (LB)` = round(`APR (LB)`, 1),
    `APR (UB)` = round(`APR (UB)`, 1),
    LPR = round(LPR, 1),
    `LPR CI LB` = round(LPR_CI_LB, 1),
    `LPR CI UB` = round(LPR_CI_UB, 1)
  )

# Save raw numeric output too, including APR CI columns.
write.csv(
  results_pct,
  file = "results/table2_raw.csv",
  row.names = FALSE
)

table2_tex <- tt(
  table2,
  caption = "Table 2"
)

print(table2_tex)

save_tt(
  table2_tex,
  file = "results/table2.tex"
)
