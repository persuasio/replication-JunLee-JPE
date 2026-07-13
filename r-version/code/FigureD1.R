############################
### R code for Figure D1 ###
############################

# Slightly modified version of the Stata replication package

rm(list = ls())

library(here)
library(readxl)
library(dplyr)
library(tidyr)

# Project directories
input_file <- here::here("output", "figureD1input.xlsx")
output_file <- here::here("output", "figureD1.pdf")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Read Figure D1 input created by figureD1input.R
figureD1_data <- read_excel(
  input_file,
  col_names = TRUE
) |>
  dplyr::select(where(is.numeric)) |>
  tidyr::drop_na()

stopifnot(ncol(figureD1_data) == 2L)

v_grid <- seq(0.40, 0.60, by = 0.01)

avg_unity <- figureD1_data[[1]][1]
avg_ovr   <- figureD1_data[[2]][1]

mpr_unity <- figureD1_data[[1]][-1]
mpr_ovr   <- figureD1_data[[2]][-1]

pdf(
  file = output_file,
  width = 10,
  height = 5
)

par(
  mfrow = c(1, 2),
  mar = c(4.5, 4.5, 3, 1)
)

plot(
  v_grid,
  mpr_unity,
  type = "l",
  lty = "solid",
  col = "blue",
  xlab = "v",
  ylab = "Marginal persuasion rate",
  main = "Not Vote for Unity"
)

abline(
  h = avg_unity,
  lty = "dashed",
  col = "black"
)

legend(
  "topright",
  legend = c("Marginal", "Average"),
  lty = c("solid", "dashed"),
  col = c("blue", "black"),
  bty = "n"
)

plot(
  v_grid,
  mpr_ovr,
  type = "l",
  lty = "solid",
  col = "blue",
  xlab = "v",
  ylab = "Marginal persuasion rate",
  main = "Vote for OVR"
)

abline(
  h = avg_ovr,
  lty = "dashed",
  col = "black"
)

legend(
  "topleft",
  legend = c("Marginal", "Average"),
  lty = c("solid", "dashed"),
  col = c("blue", "black"),
  bty = "n"
)

dev.off()

message("Saved Figure D1 to: ", output_file)
