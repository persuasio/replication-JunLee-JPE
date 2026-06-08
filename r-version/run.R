# run.R
# **********************
# OVERVIEW
#   This script generates selected tables for the paper:
#       "IDENTIFYING THE EFFECT OF PERSUASION" (Sung Jae Jun and Sokbae Lee)
#   All data are stored in /data
#   All results are outputted to /results
#
# SOFTWARE REQUIREMENTS
#   R version 4.0 or newer
#   R package: persuasio (and dependencies)
#
# TO PERFORM A CLEAN RUN, DELETE THE FOLLOWING FOLDER:
#   results/
# **********************

# **********************
# Parameters defined by user
# **********************
# Number of bootstraps for Table 2
# The paper used nbt = 10000.
# Use nbt = 10 to check whether replication code runs without an error.
nbt  <- 10
seed <- 987975

# **********************
# Setup
# **********************
library(here)        # portable paths
library(haven)       # read .dta
library(persuasio)   # your R package

set.seed(seed)

# Record session info (mirrors Stata's system parameter log)
cat("Begin:", format(Sys.time()), "\n")
cat("R version:", R.version$version.string, "\n")
cat("Platform:", R.version$platform, "\n")
sink(here("scripts/logs", paste0(format(Sys.time(), "%Y.%m.%d-%H.%M.%S"), ".log.txt")),
     split = TRUE)

# Create output directories if they don't exist
dir.create(here("results/tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(here("results/figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("scripts/logs"),    recursive = TRUE, showWarnings = FALSE)

# Load all datasets once, pass as needed
source(here("R/utils/load_data.R"))
data <- load_data()

# **********************
# Run analysis
# **********************

# TABLE 1. Persuasion Rates: Papers on Voter Turnout
source(here("R/tables/table1.R"))

# CY19_data.R already run; ChenYang2019.dta is pre-built (mirrors commented-out do file)
# source(here("R/data_prep/CY19_data.R"))

# TABLE 2. Persuasion Rates of Exposure to Uncensored Internet
source(here("R/tables/table2.R"))   # uses nbt and seed defined above

# TABLE 3. Summary Statistics of the GKB Data
source(here("R/tables/table3.R"))

# TABLE 4. Estimates of the Key Parameters
source(here("R/tables/table4.R"))

# TABLE D1. Persuasion Rates: Fox News Effects
source(here("R/tables/tableD1.R"))

# FIGURE D1. Estimates of Marginal and Average Persuasion Rates
#   figureD1input.R generates the input object; FigureD1.R renders the plot
source(here("R/figures/figureD1input.R"))
source(here("R/figures/FigureD1.R"))

# TABLE E1. Persuasive Effect by Treatment in Landry et al. (2006)
source(here("R/tables/tableE1.R"))

# TABLE E2. Persuasive Effect by Treatment in DLM
source(here("R/tables/tableE2.R"))

# TABLE H1. Persuasion Rates: NTV Effects Using a Binary Instrument
source(here("R/tables/tableH1.R"))

# **********************
cat("End:", format(Sys.time()), "\n")
sink()
