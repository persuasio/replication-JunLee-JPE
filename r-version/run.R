# **********************
# * SCRIPT:   run.R
# * PURPOSE:  Runs selected replication scripts for
# *           "IDENTIFYING THE EFFECT OF PERSUASION"
# *           (Sung Jae Jun and Sokbae Lee)
# **********************

# Stata analogue:
#   global Persuasion ".../replication-JunLee-JPE-main"
#   global nbt = 10
#   set seed 987975
#   cap mkdir "$Persuasion/scripts/logs"
#   cap mkdir "$Persuasion/results"

# -----------------------------------------------------------------------------
# User parameters
# -----------------------------------------------------------------------------

# Root directory containing /data, /scripts, and /results.
# Significant change from Stata: R uses a normal object instead of a global macro.
# If this script is run from the repository root, getwd() should be correct.
Persuasion <- getwd()

# Number of bootstraps for Table 2.
# The Stata trace uses nbt = 10 for a quick run; the paper used 10000.
nbt <- 10

# Seed for replicability.
set.seed(987975)

# -----------------------------------------------------------------------------
# Directory setup
# -----------------------------------------------------------------------------

if (missing(Persuasion) || is.na(Persuasion) || Persuasion == "") {
  stop("Persuasion root directory is not defined.")
}

scripts_dir <- file.path(Persuasion, "scripts")
results_dir <- file.path(Persuasion, "results")
logs_dir <- file.path(scripts_dir, "logs")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Some translated scripts write to outreg_dir because that naming convention was
# used in prior CY19 translations. Alias it to results_dir for consistency.
outreg_dir <- results_dir

# -----------------------------------------------------------------------------
# Package setup
# -----------------------------------------------------------------------------

required_packages <- c(
  "haven",
  "dplyr",
  "tibble",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "broom",
  "sandwich",
  "lmtest",
  "fixest",
  "tinytable",
  "persuasio"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running replication scripts: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# Source helper
# -----------------------------------------------------------------------------

run_script <- function(script_name) {

  script_path <- file.path(scripts_dir, script_name)

  if (!file.exists(script_path)) {
    warning("Skipping missing script: ", script_path)
    return(invisible(FALSE))
  }

  message("Running ", script_name, " ...")

  # Significant change from Stata: source() evaluates in the current R session,
  # roughly like Stata do-files sharing globals/macros. local = FALSE preserves
  # objects such as Persuasion, nbt, results_dir, and outreg_dir.
  source(script_path, local = FALSE)

  message("Finished ", script_name)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Script execution order
# -----------------------------------------------------------------------------

# This order follows the run.do workflow implied by the trace log and the
# translated files created from the individual trace logs. If the original
# run.do includes extra scripts not yet translated, add them here.
replication_scripts <- c(
  "table1.R",
  "table2.R",
  "table3.R",
  "table4.R",
  "tableD1.R",
  "tableD2.R",
  "figureD1input.R",
  "tableE1.R",
  "tableE2.R",
  "tableH1.R"
)

run_status <- tibble::tibble(
  script = replication_scripts,
  ran = purrr::map_lgl(replication_scripts, run_script)
)

write.csv(
  run_status,
  file.path(logs_dir, "run_status.csv"),
  row.names = FALSE
)

message("Replication run complete. Results directory: ", results_dir)
