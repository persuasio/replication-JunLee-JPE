# ************
# * SCRIPT:   table3.R
# * PURPOSE:  Creates Table 3
# *
# * ACKNOWLEDGMENT
# *       The original dataset "publicdata.dta" is from
# *       Gerber, Karlan, and Bergan (2009, AEJ Applied).
# ************

library(haven)
library(dplyr)
library(tidyr)
library(tibble)
library(tinytable)

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

# The Stata replication file also requires survey == 1.
# However, `survey` is not included in `publicdata.dta`.
# Requiring nonmissing follow-up outcome and readership variables
# reproduces the same 701-observation estimation sample.

gkb <- read_dta(data_path) %>%
  filter(
    .data$times != 1,
    !is.na(.data$voteddem_all),
    !is.na(.data$readsome)
  )

# Reproduce one Table 3 panel
make_table3 <- function(data, z_value, panel_label) {
  panel_body <- data %>%
    filter(.data$post == z_value) %>%
    count(
      .data$voteddem_all,
      .data$readsome,
      name = "n"
    ) %>%
    complete(
      voteddem_all = 0:1,
      readsome = 0:1,
      fill = list(n = 0L)
    ) %>%
    arrange(.data$voteddem_all, .data$readsome) %>%
    pivot_wider(
      names_from = .data$readsome,
      values_from = .data$n,
      names_prefix = "T="
    ) %>%
    mutate(
      Y = paste0("Y=", .data$voteddem_all),
      Total = .data$`T=0` + .data$`T=1`
    ) %>%
    select(
      Y,
      `T=0`,
      `T=1`,
      Total
    )

  panel_total <- panel_body %>%
    summarise(
      Y = "Total",
      `T=0` = sum(.data$`T=0`),
      `T=1` = sum(.data$`T=1`),
      Total = sum(.data$Total)
    )

  bind_rows(
    tibble(
      Y = panel_label,
      `T=0` = NA_integer_,
      `T=1` = NA_integer_,
      Total = NA_integer_
    ),
    panel_body,
    panel_total
  )
}

table3_data <- bind_rows(
  make_table3(
    data = gkb,
    z_value = 1,
    panel_label = "Washington Post (Z = 1):"
  ),
  make_table3(
    data = gkb,
    z_value = 0,
    panel_label = "Control (Z = 0):"
  )
)

# Save results.
write.csv(
  table3_data,
  file = file.path(
    output_dir,
    "table3_raw.csv"
  ),
  row.names = FALSE,
  na = ""
)

table3_display <- table3_data %>%
  mutate(
    across(
      c(`T=0`, `T=1`, Total),
      ~ ifelse(is.na(.), "", as.character(.))
    )
  ) %>%
  rename(
    `Voted for Democrat` = Y,
    `T_i = 0` = `T=0`,
    `T_i = 1` = `T=1`
  )

write.csv(
  table3_display,
  file = file.path(
    output_dir,
    "table3_display.csv"
  ),
  row.names = FALSE,
  na = ""
)

table3_tex <- tt(
  table3_display,
  caption = paste(
    "Summary Statistics of Data from",
    "Gerber, Karlan, and Bergan (2009)"
  )
)

print(table3_tex)

save_tt(
  table3_tex,
  output = file.path(
    output_dir,
    "table3.tex"
  ),
  overwrite = TRUE
)
