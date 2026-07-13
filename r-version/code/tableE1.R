# ************
# * SCRIPT:   tableE1.do
# * PURPOSE:  Creates Table E1
# *
# * ACKNOWLEDGMENT
# *       The summary statistic used here is from Landry, Lange, List, Price, and Rupp (QJE, 2006).
# *   Specifically, the first three rows of Table II in their paper provide information to compute
# *   the following quantities:
# *   - P(Y=1|Z=1) is obtained by (# of households that contributed)/(Total households approached);
# *   - e(1) = P(Z=1) is obtained by (Total households home)/(Total households approached).
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
# * SCRIPT:   tableE1.R
# * PURPOSE:  Creates Table E1
# ************

# Original Stata source: trace_tableE1.txt.
# This table uses published summary counts, not microdata.

landry_et_al <- matrix(
  c(
    1186, 446, 113,
    1282, 453,  67,
    963, 363, 165,
    1402, 493, 177
  ),
  ncol = 3,
  byrow = TRUE
)

# Stata: append "All" row by summing rows.
landry_et_al <- rbind(
  landry_et_al,
  colSums(landry_et_al)
)

results <- vector("list", nrow(landry_et_al))

for (j in seq_len(nrow(landry_et_al))) {

  pr_y1_z1 <- landry_et_al[j, 3] / landry_et_al[j, 1]
  pr_z1    <- landry_et_al[j, 2] / landry_et_al[j, 1]

  # Stata formulas:
  # theta_lb    = P(Y=1|Z=1)
  # theta_ub    = P(Y=1|Z=1) + 1 - P(Z=1)
  # theta_local = P(Y=1|Z=1) / P(Z=1)
  results[[j]] <- tibble(
    `P(Y=1|Z=1)` = pr_y1_z1,
    `e(1)`       = pr_z1,
    `APR (LB)`   = pr_y1_z1,
    `APR (UB)`   = pr_y1_z1 + 1 - pr_z1,
    LPR          = pr_y1_z1 / pr_z1
  )
}

tableE1_data <- bind_rows(results) * 100

tableE1_data <- bind_cols(
  Treatment = c(
    "VCM",
    "VCM with seed money",
    "Single-prize lottery",
    "Multiple-prize lottery",
    "All"
  ),
  as_tibble(tableE1_data)
)

tableE1_data[-1] <- lapply(tableE1_data[-1], round, digits = 1)

tableE1 <- tt(
  tableE1_data,
  caption = "Persuasive Effect by Treatment in Landry et al. (2006)"
)

save_tt(
  tableE1,
  file = file.path(results_dir, "tableE1.tex")
)

print(tableE1)
