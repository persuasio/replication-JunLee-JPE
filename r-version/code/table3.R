# ************
# * SCRIPT:   table3.do
# * PURPOSE:  Creates Table 3
# *
# * ACKNOWLEDGMENT
# *       The orginal dataset "publicdata.dta" is from Gerber, Karlan, and Berg an (2009, AEJ Applied).
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
# * SCRIPT:   table3.R
# * PURPOSE:  Creates Table 3
# ************

# Original Stata source: trace_table3.txt.
# Main workflow:
#   1. load GKB publicdata
#   2. drop Washington Times treatment
#   3. keep estimation sample
#   4. tabulate Y by T separately for Z = 1 and Z = 0

gkb <- read_dta(
  file.path(
    persuasion_dir,
    "data",
    "GerberKarlanBergan2009",
    "publicdata.dta"
  )
)

gkb <- gkb %>%
  filter(times != 1) %>%
  mutate(
    data_avail = as.integer(
      survey == 1 &
        !is.na(voteddem_all) &
        !is.na(readsome)
    )
  ) %>%
  filter(data_avail == 1)

# Stata: tab voteddem_all readsome if post == 1, matcell(treat)
treat <- gkb %>%
  filter(post == 1) %>%
  count(voteddem_all, readsome, name = "n") %>%
  tidyr::complete(
    voteddem_all = 0:1,
    readsome = 0:1,
    fill = list(n = 0)
  ) %>%
  arrange(voteddem_all, readsome) %>%
  tidyr::pivot_wider(
    names_from = readsome,
    values_from = n,
    names_prefix = "T="
  ) %>%
  mutate(Y = paste0("Y=", voteddem_all)) %>%
  select(Y, `T=0`, `T=1`)

table3a <- tt(
  treat,
  caption = "Table 3: The Washington Post (Z = 1)"
)

save_tt(
  table3a,
  file = file.path(results_dir, "table3a.tex")
)

# Stata: tab voteddem_all readsome if post == 0, matcell(control)
control <- gkb %>%
  filter(post == 0) %>%
  count(voteddem_all, readsome, name = "n") %>%
  tidyr::complete(
    voteddem_all = 0:1,
    readsome = 0:1,
    fill = list(n = 0)
  ) %>%
  arrange(voteddem_all, readsome) %>%
  tidyr::pivot_wider(
    names_from = readsome,
    values_from = n,
    names_prefix = "T="
  ) %>%
  mutate(Y = paste0("Y=", voteddem_all)) %>%
  select(Y, `T=0`, `T=1`)

table3b <- tt(
  control,
  caption = "Table 3: Control (Z = 0)"
)

save_tt(
  table3b,
  file = file.path(results_dir, "table3b.tex")
)

print(table3a)
print(table3b)
