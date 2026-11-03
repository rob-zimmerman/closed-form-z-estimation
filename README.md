# Code for *An Asymptotic Framework for Closed-Form Estimators from Surrogate Estimating Equations*

This repository contains the simulation and application code accompanying the
paper. It reproduces the twelve figures for the gamma, Zeta, von Mises, and
generalized Pareto examples, the black cherry tree analysis, and the
high-precision numerical check for the gamma likelihood equations.

## Requirements

- R 4.2 or later
- R packages `dplyr`, `ggplot2`, `scales`, and `tidyr`
- Ghostscript (`gs`) for font-embedded PDF figures
- Python 3 for the gamma numerical audit
- CMU Serif for figures matching the manuscript typography

Install the R packages with:

```r
install.packages(c("dplyr", "ggplot2", "scales", "tidyr"))
```

If CMU Serif is unavailable, set the `FIGURE_FONT` environment variable to an
installed font before running the scripts.

## Running the analyses

Run all analyses from the repository root:

```sh
Rscript scripts/run_all.R
```

The manuscript run uses 2,000 Monte Carlo replications at every design point
and may take substantial time. A short smoke test exercises every estimator and
figure on a reduced grid:

```sh
Rscript scripts/run_all.R --quick
```

Figures are written to `output/figures/` and compact numerical summaries to
`output/summaries/`. Generated output is ignored by Git.

## Files

- `scripts/gamma_simulation.R`: known-rate gamma, unknown-rate gamma, and gamma
  regression simulations; Figures 1--3.
- `scripts/gamma_trees.R`: black cherry tree gamma-regression application.
- `scripts/zeta_simulation.R`: Zeta simulations; Figures 4--6.
- `scripts/von_mises_simulation.R`: von Mises simulations; Figures 7--9.
- `scripts/gpd_simulation.R`: generalized Pareto simulations; Figures 10--12.
- `scripts/verify_gamma_numerics.py`: 90-digit audit of the large-shape gamma
  roots.

Every simulation cell uses a fixed cell-specific seed. The scripts create only
the figures in the paper and the summary tables needed to reproduce its
reported numerical values; they do not save individual Monte Carlo
replications.
