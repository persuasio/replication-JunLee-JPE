# ************
# * SCRIPT:   figureD1input.R
# * PURPOSE:  Creates input for Figure D1
# *
# * ACKNOWLEDGMENT
# *       The original dataset "NTV_Individual_Data.dta" is from
# *       Enikolopov, Petrova, and Zhuravskaya (AER, 2011).
# ************

library(haven)
library(dplyr)
library(tibble)
library(tidyr)
library(tinytable)
library(writexl)
library(here)

output_dir <- here::here("output")
data_path <- here::here(
  "data",
  "EnikolopovPetrovaZhuravskaya2011",
  "NTV_Individual_Data.dta"
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

weighted_mean_safe <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0

  if (!any(keep)) {
    stop("No observations with nonmissing values and positive weights.", call. = FALSE)
  }

  weighted.mean(x[keep], w[keep])
}

ntv <- read_dta(data_path) %>%
  mutate(
    # Reproduce the Stata initialization and replacement commands
    # Missing source outcomes remain coded as zero
    y_vote_Unity = as.integer(!is.na(vote_Unity) & vote_Unity == 0),
    y_vote_OVR = as.integer(!is.na(vote_OVR) & vote_OVR == 1)
  )

# Predict the exposure rate e(X, Z) = P(T = 1 | X, Z).
# The specification is linear in X and cubic in Z, with selected X-by-Z
# interactions, as in the Stata
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
  weights = kishweig,
  na.action = na.exclude
)

# Stata's predict leaves observations outside the estimation sample missing
ntv$phat <- clip01(predict(exposure_fit, newdata = ntv))

mte_grid <- seq(0.40, 0.60, by = 0.01)

# Reproduce margins, dydx(phat) at(phat = ...): calculate each observation's
# derivative and then average it using kishweig over the model sample
marginal_effect_at_phat <- function(fit, data, grid = mte_grid) {
  b <- coef(fit)
  estimation_sample <- rownames(data) %in% rownames(model.frame(fit))

  vapply(grid, function(v) {
    derivative <- rep(0, nrow(data))

    derivative <- derivative + ifelse(
      "phat" %in% names(b),
      b[["phat"]],
      0
    )
    derivative <- derivative + ifelse(
      "I(phat^2)" %in% names(b),
      2 * b[["I(phat^2)"]] * v,
      0
    )
    derivative <- derivative + ifelse(
      "I(phat^3)" %in% names(b),
      3 * b[["I(phat^3)"]] * v^2,
      0
    )

    for (x in c("male", "age", "educ1", "married", "consump")) {
      term_names <- c(
        paste0("phat:", x),
        paste0(x, ":phat")
      )
      matched_term <- term_names[term_names %in% names(b)]
      interaction_coefficient <- if (length(matched_term) > 0L) {
        b[[matched_term[[1L]]]]
      } else {
        0
      }

      derivative <- derivative + interaction_coefficient * data[[x]]
    }

    derivative[!estimation_sample] <- NA_real_
    weighted_mean_safe(derivative, data$kishweig)
  }, numeric(1))
}

make_figureD1input_column <- function(data, party) {
  y_var <- paste0("y_vote_", party)
  notwatch_var <- paste0("notwatch_vote_", party)

  data[[notwatch_var]] <- as.integer(
    data[[y_var]] == 1 &
      !is.na(data$Watches_NTV_1999) &
      data$Watches_NTV_1999 == 0
  )

  rhs <- paste(
    c(
      "male", "age", "educ1", "married", "consump",
      "logpop98", "wage98_ln",
      "phat", "I(phat^2)", "I(phat^3)",
      "phat:male", "phat:age", "phat:educ1",
      "phat:married", "phat:consump"
    ),
    collapse = " + "
  )

  y_fit <- lm(
    as.formula(paste(y_var, "~", rhs)),
    data = data,
    weights = kishweig,
    na.action = na.exclude
  )

  notwatch_fit <- lm(
    as.formula(paste(notwatch_var, "~", rhs)),
    data = data,
    weights = kishweig,
    na.action = na.exclude
  )

  num <- marginal_effect_at_phat(y_fit, data)
  den <- 1 + marginal_effect_at_phat(notwatch_fit, data)
  mpr <- num / den

  # Stata computes rowsum(num) / rowsum(den), not mean(num / den).
  average_mpr <- sum(num) / sum(den)

  tibble(
    row = c(
      "Avg between 0.4 and 0.6",
      sprintf("v = %.2f", mte_grid)
    ),
    value = c(average_mpr, mpr)
  )
}

unity_results <- make_figureD1input_column(ntv, "Unity")
ovr_results <- make_figureD1input_column(ntv, "OVR")

# Raw numeric matrix in the same row/column orientation as the Stata table
figureD1input_matrix <- unity_results %>%
  select(row, unity = value) %>%
  left_join(
    ovr_results %>%
      select(row, ovr = value),
    by = "row"
  ) %>%
  transmute(
    Estimand = row,
    `Not Vote for Unity` = unity,
    `Vote for OVR` = ovr
  )

# Display version: Stata reports three decimal places.
figureD1input_data <- figureD1input_matrix %>%
  mutate(
    across(
      -Estimand,
      ~ sprintf("%.3f", .x)
    )
  )


# Save results using the same output convention as tableD1.R
write.csv(
  figureD1input_matrix,
  file = file.path(
    output_dir,
    "figureD1input_raw.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  figureD1input_data,
  file = file.path(
    output_dir,
    "figureD1input_display.csv"
  ),
  row.names = FALSE,
  na = ""
)

writexl::write_xlsx(
  figureD1input_matrix,
  path = file.path(
    output_dir,
    "figureD1input.xlsx"
  )
)

figureD1input_tex <- tt(
  figureD1input_data,
  caption = "Estimates of Marginal and Average Persuasion Rates"
)

print(figureD1input_tex)

save_tt(
  figureD1input_tex,
  output = file.path(
    output_dir,
    "figureD1input.tex"
  ),
  overwrite = TRUE
)
