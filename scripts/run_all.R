#!/usr/bin/env Rscript

# Run every simulation and the black cherry tree application from the
# repository root. Pass --quick to use a small grid for a short smoke test.

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Unable to determine the location of scripts/run_all.R")
}
script_path <- normalizePath(sub("^--file=", "", file_argument))
repository_root <- dirname(dirname(script_path))
setwd(repository_root)

if (quick) {
  # Small replication counts test every estimator and figure without attempting
  # to reproduce the Monte Carlo precision of the manuscript.
  Sys.setenv(
    GAMMA_UNIFIED_QUICK = "1",
    GAMMA_UNIFIED_B = "40",
    GAMMA_UNIFIED_CHUNK_SIZE = "20",
    ZETA_SIM_QUICK = "1",
    ZETA_SIM_B = "40",
    ZETA_SIM_CHUNK = "20",
    VONMISES_SIM_QUICK = "1",
    VONMISES_SIM_B = "40",
    VONMISES_SIM_CHUNK = "20",
    GPD_SIM_QUICK = "1",
    GPD_SIM_B = "30",
    GPD_SIM_CHUNK = "15"
  )
}

run_command <- function(command, arguments, label) {
  message("\n== ", label, " ==")
  status <- system2(command, arguments)
  if (!identical(status, 0L)) {
    stop(label, " failed with exit status ", status)
  }
}

rscript <- file.path(R.home("bin"), "Rscript")
simulation_scripts <- c(
  "Gamma simulation" = "scripts/gamma_simulation.R",
  "Zeta simulation" = "scripts/zeta_simulation.R",
  "Von Mises simulation" = "scripts/von_mises_simulation.R",
  "GPD simulation" = "scripts/gpd_simulation.R"
)

for (label in names(simulation_scripts)) {
  run_command(rscript, simulation_scripts[[label]], label)
}

run_command(rscript, "scripts/gamma_trees.R", "Black cherry tree application")

python <- Sys.which("python3")
if (!nzchar(python)) {
  stop("python3 is required for the high-precision gamma numerical audit")
}
run_command(
  python,
  "scripts/verify_gamma_numerics.py",
  "High-precision gamma numerical audit"
)

message("\nAll analyses completed successfully.")
