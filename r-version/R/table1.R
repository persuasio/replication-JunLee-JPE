# ************
# * SCRIPT:   table1.R
# * PURPOSE:  Creates Table 1
# ************

library(tibble)
library(persuasio)
library(tinytable)

# Input for Table 1
turnout <- matrix(
  c(
    0.472, 0.448, 0.279, 0,
    0.310, 0.286, 0.293, 0,
    0.711, 0.660, 0.737, 0,
    0.416, 0.405, 0.414, 0,
    0.455, 0.435, 0.800, 0,
    0.700, 0.690, 0.250, 0
  ),
  ncol = 4,
  byrow = TRUE
)

reference <- c(
  "Green and Gerber (2000)",
  "Green, Gerber, and Nickerson (2003)",
  "Green and Gerber (2001)",
  "Green and Gerber (2001)",
  "Gentzkow (2006)",
  "Gentzkow, Shapiro, and Sinkinson (2011)"
)

results <- vector("list", nrow(turnout))
for (j in seq_len(nrow(turnout))) {
  y1 <- turnout[j, 1]
  y0 <- turnout[j, 2]
  e1 <- turnout[j, 3]
  e0 <- turnout[j, 4]

  out <- calc4persuasio(
    y1 = y1,
    y0 = y0,
    e1 = e1,
    e0 = e0
  )

  persuasion_DK <- (y1 - y0) / ((e1 - e0) * (1 - y0))

  results[[j]] <- c(
    DK     = persuasion_DK,
    APR_LB = out$apr[1],
    APR_UB = out$apr[2],
    LPR_LB = out$lpr[1],
    LPR_UB = out$lpr[2]
  )
}

results <- as.data.frame(do.call(rbind, results))

table1 <- tibble(
  Study = reference,
  DK = results$DK,
  `APR (LB)` = results$APR_LB,
  `APR (UB)` = results$APR_UB,
  `LPR (LB)` = results$LPR_LB,
  `LPR (UB)` = results$LPR_UB
)

table1[-1] <- lapply(table1[-1], round, digits = 3)

table1 <- tt(
  table1,
  caption = "Persuasion Rates: Papers on Voter Turnout"
)
print(table1)

save_tt(
  table1,
  file = "results/table1.tex"
)
