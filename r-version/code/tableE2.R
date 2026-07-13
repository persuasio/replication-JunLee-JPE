# ************
# * SCRIPT:   tableE2.R
# * PURPOSE:  Creates Table E2
# *
# * ACKNOWLEDGMENT
# *       The original dataset "CharityOutputQJE" is from
# *       DellaVigna, List, and Malmendier (QJE, 2012).
# ************

library(haven)
library(dplyr)
library(tibble)
library(tinytable)
library(here)

output_dir <- here::here("output")
data_path <- here::here(
  "data",
  "DellaVignaListMalmendier2012",
  "CharityOutputQJE.dta"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ",
    data_path,
    call. = FALSE
  )
}

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

common_means <- function(data, variables) {
  keep <- complete.cases(data[, variables, drop = FALSE])

  if (!any(keep)) {
    stop(
      "No common nonmissing observations for: ",
      paste(variables, collapse = ", "),
      call. = FALSE
    )
  }

  vapply(
    data[keep, variables, drop = FALSE],
    mean,
    numeric(1)
  )
}

charity <- read_dta(data_path) %>%
  mutate(
    treatment = if_else(treatment == "2Ww", "W", treatment)
  ) %>%
  filter(nosolsign != 1) %>%
  select(-nosolsign) %>%
  filter(toeliminate != 1) %>%
  select(-toeliminate) %>%
  filter(
    !(solicitor == "Angelena" & date %in% as.Date(c("2008-07-27", "2008-07-13"))),
    !(solicitor == "Shedora" & date %in% as.Date(c("2008-07-27", "2008-08-10"))),
    !(solicitor == "Tehmur" & date == as.Date("2008-06-01")),
    !(solicitor == "Phillip" & date == as.Date("2008-08-09")),
    !(solicitor == "Robert" & date == as.Date("2008-07-13") & hour %in% c(11, 13))
  ) %>%
  mutate(
    month = as.integer(format(date, "%m")),
    dwave = case_when(
      year == 2008 & month %in% c(7, 8) & charity == "LaRabida" ~ 1,
      year == 2008 &
        (month %in% c(4, 5, 6) | (charity == "Ecu" & month == 7)) ~ 0,
      year == 2008 & month %in% c(9, 10) ~ 2,
      TRUE ~ NA_real_
    ),
    # Stata: egen sodate = concat(solicitor date)
    sodate = interaction(solicitor, date, drop = TRUE),
    treatmentby = interaction(treatment, charity, drop = TRUE),
    grsol = factor(solicitor),
    grdatloc = interaction(date, location, drop = TRUE),
    grdatlocsol = interaction(date, location, solicitor, drop = TRUE),
    grhour = factor(hour),
    grarea = factor(area_rank),
    amt_donate = if_else(is.na(amt_donate), 0, amt_donate),
    dW = as.integer(treatment == "W"),
    dOo = as.integer(treatment == "Oo"),
    dWEcu = as.integer(treatment == "W" & charity == "Ecu"),
    dWLar = as.integer(treatment == "W" & charity == "LaRabida"),
    dOoEcu = as.integer(treatment == "Oo" & charity == "Ecu"),
    dOoLar = as.integer(treatment == "Oo" & charity == "LaRabida"),
    dEcu = as.integer(charity == "Ecu")
  )

# Reproduce Stata's treatment-name cleanup.
for (x in c("0d5m", "0d10m", "5d5m", "5d10m", "10d10m", "10d5m")) {
  for (y in c("Nw", "W", "Oo")) {
    old_value <- paste0(y, "-", x)
    new_value <- paste0(y, x)
    charity$treatment[charity$treatment == old_value] <- new_value
  }
}

# These indicators are generated in the Stata source before the final models.
for (x in c("Nw0d10m", "W0d10m", "W0d5m", "W10d10m")) {
  charity[[paste0("d", x, "08")]] <- as.integer(
    charity$treatment == x & charity$year == 2008
  )
}

for (x in c(
  "Nw0d5m", "Nw5d5m", "W0d10m", "W0d5m",
  "W10d5m", "W5d5m", "Oo0d5m", "Oo5d5m"
)) {
  charity[[paste0("d", x, "09")]] <- as.integer(
    charity$treatment == x & charity$year == 2009
  )
}

# Stata drops the omitted categories.
charity <- charity %>%
  select(-any_of(c("dNw0d10m08", "dNw0d5m09")))

make_tableE2_rows <- function(data, charity_name) {
  tmp <- data %>%
    filter(charity == charity_name)

  if (nrow(tmp) == 0L) {
    stop("No observations found for charity: ", charity_name, call. = FALSE)
  }

  # Stata: xi: reg answer dW dOo i.grsol i.grdatloc i.grhour i.grarea
  a_fit <- lm(
    answer ~ dW + dOo + grsol + grdatloc + grhour + grarea,
    data = tmp,
    na.action = na.exclude
  )

  coeff_a_W <- unname(coef(a_fit)[["dW"]])
  coeff_a_Oo <- unname(coef(a_fit)[["dOo"]])
  ahat <- predict(a_fit, newdata = tmp)

  tmp <- tmp %>%
    mutate(
      ahat_Nw = clip01(ahat - coeff_a_W * dW - coeff_a_Oo * dOo),
      ahat_W = clip01(ahat_Nw + coeff_a_W),
      ahat_Oo = clip01(ahat_Nw + coeff_a_Oo)
    )

  # Stata: xi: reg saidyes dW dOo i.grsol i.grdatloc i.grhour i.grarea
  y_fit <- lm(
    saidyes ~ dW + dOo + grsol + grdatloc + grhour + grarea,
    data = tmp,
    na.action = na.exclude
  )

  coeff_y_W <- unname(coef(y_fit)[["dW"]])
  coeff_y_Oo <- unname(coef(y_fit)[["dOo"]])
  yhat <- predict(y_fit, newdata = tmp)

  tmp <- tmp %>%
    mutate(
      yhat_Nw = clip01(yhat - coeff_y_W * dW - coeff_y_Oo * dOo),
      yhat_W = clip01(yhat_Nw + coeff_y_W),
      yhat_Oo = clip01(yhat_Nw + coeff_y_Oo),
      lb_Nw = yhat_Nw,
      lb_W = yhat_W,
      lb_Oo = yhat_Oo,
      ub_Nw = yhat_Nw + 1 - ahat_Nw,
      ub_W = yhat_W + 1 - ahat_W,
      ub_Oo = yhat_Oo + 1 - ahat_Oo
    )

  arm_labels <- c(
    Nw = "Baseline",
    W = "Flyer",
    Oo = "Opt-Out"
  )

  bind_rows(lapply(names(arm_labels), function(arm) {
    variables <- c(
      paste0("yhat_", arm),
      paste0("ahat_", arm),
      paste0("lb_", arm),
      paste0("ub_", arm)
    )

    means <- common_means(tmp, variables)

    tibble(
      Row = paste(
        ifelse(charity_name == "LaRabida", "La Rabida", charity_name),
        arm_labels[[arm]],
        sep = ": "
      ),
      `P(Y=1|Z=1)` = unname(means[[1]]),
      `e(1)` = unname(means[[2]]),
      `APR (LB)` = unname(means[[3]]),
      `APR (UB)` = unname(means[[4]]),
      LPR = unname(means[[1]] / means[[2]])
    )
  }))
}

# Stata stacks La Rabida first and Ecu second.
tableE2_matrix <- bind_rows(
  make_tableE2_rows(charity, "LaRabida"),
  make_tableE2_rows(charity, "Ecu")
)

# Stata multiplies the complete result matrix by 100.
tableE2_raw <- tableE2_matrix %>%
  mutate(across(-Row, ~ 100 * .x))

# Display version: Stata uses one decimal place.
tableE2_data <- tableE2_raw %>%
  mutate(across(-Row, ~ sprintf("%.1f", .x)))

stopifnot(
  nrow(tableE2_raw) == 6L,
  identical(
    tableE2_raw$Row,
    c(
      "La Rabida: Baseline",
      "La Rabida: Flyer",
      "La Rabida: Opt-Out",
      "Ecu: Baseline",
      "Ecu: Flyer",
      "Ecu: Opt-Out"
    )
  )
)

write.csv(
  tableE2_raw,
  file = file.path(output_dir, "tableE2_raw.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  tableE2_data,
  file = file.path(output_dir, "tableE2_display.csv"),
  row.names = FALSE,
  na = ""
)

tableE2_tex <- tt(
  tableE2_data,
  caption = "Persuasive Effect by Treatment in DLM"
)

print(tableE2_tex)

save_tt(
  tableE2_tex,
  output = file.path(output_dir, "tableE2.tex"),
  overwrite = TRUE
)
