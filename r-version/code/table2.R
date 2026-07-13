# ************
# * SCRIPT:   table2.R
# * PURPOSE:  Creates Table 2
# ************

library(dplyr)
library(purrr)
library(tibble)
library(persuasio)
library(tinytable)
library(here)

output_dir <- here::here("output")
data_path <- here::here("output", "ChenYang2019.rds")

if (!file.exists(data_path)) {
  message("Output file not found. Automatically sourcing CY19_data.R to build it...")
  source(here::here("code", "CY19_data.R"))
}

cy <- readRDS(data_path)

# Use 10,000 for the final replication.
nbt <- 50

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Estimate the complier share

active_fit <- lm(
  active_user ~
    treatment_vpnonly +
    treatment_nlonly +
    treatment_vpnnl,
  data = cy
)

share_compliers <- unname(
  coef(active_fit)[["treatment_vpnnl"]]
)

# Helper: produce the Stata-style rows for one outcome

make_table2 <- function(data, y, y0 = NULL) {

  rf_data <- data %>%
    filter(
      !is.na(.data[[y]]),
      !is.na(treatment_vpnonly),
      !is.na(treatment_nlonly),
      !is.na(treatment_vpnnl)
    )

  rf_fit <- lm(
    reformulate(
      c(
        "treatment_vpnonly",
        "treatment_nlonly",
        "treatment_vpnnl"
      ),
      response = y
    ),
    data = rf_data
  )

  b_vpnnl <- unname(coef(rf_fit)[["treatment_vpnnl"]])

  if (!is.null(y0)) {
    baseline_share <- data %>%
      filter(
        treatment_vpnnl == 1,
        !is.na(.data[[y0]])
      ) %>%
      summarise(value = mean(.data[[y0]])) %>%
      pull(value)
  } else {
    baseline_share <- unname(coef(rf_fit)[["(Intercept)"]])
  }

  persuasion_DK <- b_vpnnl /
    (share_compliers * (1 - baseline_share))

  pers_data <- data %>%
    filter(treatment_main %in% c(1, 3)) %>%
    filter(
      if_all(
        all_of(c(y, "active_user", "treatment_vpnnl")),
        ~ !is.na(.x)
      )
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

  bind_rows(
    tibble(
      variable = y,
      row_type = "estimate",
      CY19 = persuasion_DK,
      APR_LB = apr_out$lb_coef,
      APR_UB = apr_out$ub_coef,
      LPR_LB = lpr_out$lpr,
      LPR_UB = NA_real_
    ),
    tibble(
      variable = y,
      row_type = "confidence interval",
      CY19 = NA_real_,
      APR_LB = apr_out$ci_lb,
      APR_UB = apr_out$ci_ub,
      LPR_LB = lpr_out$ci_lb,
      LPR_UB = lpr_out$ci_ub
    )
  )
}

# Variables in their original Stata-loop order

panel_mappings <- c(
  "info_foreign_website_w3_p" =
    "info_foreign_website_w1_p",

  "info_freq_website_for_w3_p" =
    "info_freq_website_for_w1_p",

  "az_belief_media_value_w3_p" =
    "az_belief_media_value_w1_p",

  "az_belief_media_trust_w3_p" =
    "az_belief_media_trust_w1_p",

  "bias_domestic_w3_p" =
    "bias_domestic_w1_p",

  "bias_foreign_w3_p" =
    "bias_foreign_w1_p",

  "az_belief_media_justif_w3_p" =
    "az_belief_media_justif_w1_p",

  "bias_dom_govt_policy_t1_w3_p" =
    "bias_dom_govt_policy_t1_w1_p",

  "bias_for_govt_policy_t1_w3_p" =
    "bias_for_govt_policy_t1_w1_p"
)

nonpanel_vars <- c(
  "vpn_purchase_wmt_record_p",
  "vpn_purchase_yes_p"
)

panel_results <- imap_dfr(
  panel_mappings,
  \(y0_var, y_var) {
    make_table2(
      data = cy,
      y = y_var,
      y0 = y0_var
    )
  }
)

nonpanel_results <- map_dfr(
  nonpanel_vars,
  \(y_var) {
    make_table2(
      data = cy,
      y = y_var,
      y0 = NULL
    )
  }
)

# Table formatting
#   first 2 panel outcomes = rows 1–4
#   2 purchase outcomes    = rows 19–22
#   remaining panel rows   = rows 5–18

final_outcome_order <- c(
  "info_foreign_website_w3_p",
  "info_freq_website_for_w3_p",
  "vpn_purchase_wmt_record_p",
  "vpn_purchase_yes_p",
  "az_belief_media_value_w3_p",
  "az_belief_media_trust_w3_p",
  "bias_domestic_w3_p",
  "bias_foreign_w3_p",
  "az_belief_media_justif_w3_p",
  "bias_dom_govt_policy_t1_w3_p",
  "bias_for_govt_policy_t1_w3_p"
)

results <- bind_rows(
  panel_results,
  nonpanel_results
) %>%
  mutate(
    variable = factor(
      variable,
      levels = final_outcome_order
    ),
    row_type = factor(
      row_type,
      levels = c(
        "estimate",
        "confidence interval"
      )
    )
  ) %>%
  arrange(variable, row_type) %>%
  mutate(
    variable = as.character(variable)
  )

stopifnot(nrow(results) == 22L)

# Add human-readable outcome labels
outcome_labels <- c(
  "info_foreign_website_w3_p" =
    "A.1.2. Ranked high: foreign websites",

  "info_freq_website_for_w3_p" =
    "A.1.6. Frequency of visiting foreign websites for information",

  "vpn_purchase_wmt_record_p" =
    "A.2.1. Purchase discounted tool we offered",

  "vpn_purchase_yes_p" =
    "A.2.2. Purchase any tool",

  "az_belief_media_value_w3_p" =
    "A.3. Valuation of access to foreign media outlets",

  "az_belief_media_trust_w3_p" =
    "A.4. Trust in nondomestic media outlets",

  "bias_domestic_w3_p" =
    "A.5.1. Degree of censorship on domestic news outlets",

  "bias_foreign_w3_p" =
    "A.5.2. Degree of censorship on foreign news outlets",

  "az_belief_media_justif_w3_p" =
    "A.6. Censorship unjustified",

  "bias_dom_govt_policy_t1_w3_p" =
    "A.7.1. Domestic censorship driven by government policies",

  "bias_for_govt_policy_t1_w3_p" =
    "A.7.2. Foreign censorship driven by government policies"
)


table2_matrix <- results %>%
  transmute(
    CY19 = round(100 * CY19, 1),
    APR_LB = round(100 * APR_LB, 1),
    APR_UB = round(100 * APR_UB, 1),
    LPR_LB = round(100 * LPR_LB, 1),
    LPR_UB = round(100 * LPR_UB, 1)
  )

# Show the outcome name only on the estimate row.
table2_data <- results %>%
  mutate(
    Outcome = if_else(
      row_type == "estimate",
      unname(outcome_labels[variable]),
      ""
    )
  ) %>%
  transmute(
    Outcome,

    CY19 = if_else(
      is.na(CY19),
      "",
      sprintf("%.1f", 100 * CY19)
    ),

    `APR (LB)` = if_else(
      is.na(APR_LB),
      "",
      sprintf("%.1f", 100 * APR_LB)
    ),

    `APR (UB)` = if_else(
      is.na(APR_UB),
      "",
      sprintf("%.1f", 100 * APR_UB)
    ),

    LPR = if_else(
      is.na(LPR_LB),
      "",
      sprintf("%.1f", 100 * LPR_LB)
    ),

    ` ` = if_else(
      is.na(LPR_UB),
      "",
      sprintf("%.1f", 100 * LPR_UB)
    )
  )

# Save the results

write.csv(
  table2_matrix,
  file = file.path(
    output_dir,
    "table2_raw.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  table2_data,
  file = file.path(
    output_dir,
    "table2_display.csv"
  ),
  row.names = FALSE,
  na = ""
)

table2_tex <- tt(
  table2_data,
  caption = "Table 2"
)

print(table2_tex)

save_tt(
  table2_tex,
  output = file.path(
    output_dir,
    "table2.tex"
  ),
  overwrite = TRUE
)
