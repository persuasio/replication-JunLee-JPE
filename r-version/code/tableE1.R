# ************
# * SCRIPT:   tableE1.R
# * PURPOSE:  Creates Table E1
# *
# * ACKNOWLEDGMENT
# *       The summary statistics used here are from Landry, Lange, List,
# *       Price, and Rupp (QJE, 2006).
# *
# * NOTES
# *       The first three rows of Table II in the original paper provide
# *       the counts used to calculate:
# *       - P(Y=1|Z=1) = households contributing / households approached;
# *       - e(1)       = households home / households approached.
# ************

library(dplyr)
library(tibble)
library(tinytable)
library(here)

output_dir <- here::here("output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Published summary counts in the same order as the Stata matrix:
# households approached, households home, households contributing.
landry_counts <- tribble(
  ~Treatment,                    ~approached, ~home, ~contributed,
  "VCM",                               1186,   446,         113,
  "VCM with seed money",               1282,   453,          67,
  "Single-prize lottery",               963,   363,         165,
  "Multiple-prize lottery",            1402,   493,         177
)

# Stata appends an "All" row by summing the four treatment rows first,
# and then computes all ratios from those pooled counts.
landry_counts <- bind_rows(
  landry_counts,
  landry_counts %>%
    summarise(
      Treatment = "All",
      across(c(approached, home, contributed), sum)
    )
)

# Reproduce the five quantities in the Stata trace.
tableE1_matrix <- landry_counts %>%
  transmute(
    Treatment,
    `P(Y=1|Z=1)` = 100 * contributed / approached,
    `e(1)` = 100 * home / approached,
    `APR (LB)` = 100 * contributed / approached,
    `APR (UB)` = 100 * (
      contributed / approached + 1 - home / approached
    ),
    LPR = 100 * contributed / home
  )

# Display version: Stata multiplies by 100 and reports one decimal place.
tableE1_data <- tableE1_matrix %>%
  mutate(
    across(
      -Treatment,
      ~ sprintf("%.1f", .x)
    )
  )

stopifnot(
  nrow(tableE1_matrix) == 5L,
  identical(
    tableE1_matrix$Treatment,
    c(
      "VCM",
      "VCM with seed money",
      "Single-prize lottery",
      "Multiple-prize lottery",
      "All"
    )
  ),
  identical(
    unname(round(as.numeric(tableE1_matrix[5, -1]), 1)),
    c(10.8, 36.3, 10.8, 74.5, 29.7)
  )
)

# Save the raw numeric and formatted display data using the same convention
# as tableD1.R.
write.csv(
  tableE1_matrix,
  file = file.path(
    output_dir,
    "tableE1_raw.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  tableE1_data,
  file = file.path(
    output_dir,
    "tableE1_display.csv"
  ),
  row.names = FALSE,
  na = ""
)

tableE1_tex <- tt(
  tableE1_data,
  caption = "Table E1. Persuasive Effect by Treatment in Landry et al. (2006)"
)

print(tableE1_tex)

save_tt(
  tableE1_tex,
  output = file.path(
    output_dir,
    "tableE1.tex"
  ),
  overwrite = TRUE
)
