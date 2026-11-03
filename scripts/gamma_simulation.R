#!/usr/bin/env Rscript

# Monte Carlo experiments for the gamma models in the paper.
#
# The script covers the known-rate model, the unknown-rate model, and gamma
# regression with a log link. It produces the three gamma figures and compact
# summary tables containing the quantities reported in the manuscript.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

options(scipen = 4)

# Output locations can be overridden without changing the script.
out_dir <- Sys.getenv(
  "GAMMA_UNIFIED_OUT_DIR",
  "output/summaries/gamma"
)
fig_dir <- Sys.getenv(
  "GAMMA_UNIFIED_FIG_DIR",
  "output/figures/gamma"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Full runs use the manuscript grid and 2,000 replications per cell. Quick mode
# retains all three asymptotic regimes on a two-point sample-size grid.
B <- as.integer(Sys.getenv("GAMMA_UNIFIED_B", "2000"))
chunk_size <- as.integer(Sys.getenv("GAMMA_UNIFIED_CHUNK_SIZE", "100"))
quick <- identical(Sys.getenv("GAMMA_UNIFIED_QUICK", "0"), "1")
n_grid <- if (quick) {
  c(250L, 1000L)
} else {
  c(250L, 500L, 1000L, 2500L, 5000L, 10000L, 25000L)
}
growth_grid <- c(0.10, 0.40, 0.70)
beta_rate <- 2
beta_reg <- c(0.5, -0.75)
diagnostic_n <- max(n_grid)
diagnostic_growth <- 0.40
figure_family <- Sys.getenv("FIGURE_FONT", "CMU Serif")

blue <- "#0072B2"
orange <- "#D55E00"
green <- "#009E73"
purple <- "#CC79A7"
grey <- "#4D4D4D"
growth_colours <- c(
  "0.1" = "#0072B2",
  "0.4" = "#E69F00",
  "0.7" = "#CC79A7"
)

paper_theme <- theme_bw(base_size = 11, base_family = figure_family) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.border = element_rect(colour = "grey35", linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    strip.background = element_rect(fill = "grey94", colour = "grey55", linewidth = 0.35),
    strip.text = element_text(face = "plain", size = 10.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -2, r = 0, b = 0, l = 0),
    plot.margin = margin(5, 6, 3, 5)
  )

power_of_ten_labels <- function(x) {
  exponent <- round(log10(x))
  parse(text = ifelse(exponent == 0, "1", paste0("10^", exponent)))
}

# Direct special-function evaluations lose relative accuracy when the shape is
# large. These Bernoulli expansions evaluate the same score and derivative
# without subtracting nearly equal quantities.
profile_score_direct <- function(a) {
  log(a) - digamma(a)
}

profile_score_prime_direct <- function(a) {
  1 / a - trigamma(a)
}

profile_score <- function(a) {
  answer <- profile_score_direct(a)
  large <- is.finite(a) & a >= 50
  if (any(large)) {
    x <- 1 / a[large]
    answer[large] <-
      x / 2 + x^2 / 12 - x^4 / 120 + x^6 / 252 -
      x^8 / 240 + x^10 / 132 - 691 * x^12 / 32760 + x^14 / 12
  }
  answer
}

profile_score_prime <- function(a) {
  answer <- profile_score_prime_direct(a)
  large <- is.finite(a) & a >= 50
  if (any(large)) {
    x <- 1 / a[large]
    answer[large] <-
      -x^2 / 2 - x^3 / 6 + x^5 / 30 - x^7 / 42 +
      x^9 / 30 - 5 * x^11 / 66 + 691 * x^13 / 2730 - 7 * x^15 / 6
  }
  answer
}

# Safeguarded Newton solvers for the scalar likelihood equations. Invalid
# proposals are replaced by a positive half-step before the next iteration.
inv_digamma <- function(y) {
  a <- ifelse(y >= -2.22, exp(y) + 0.5, -1 / (y + digamma(1)))
  a <- pmax(a, .Machine$double.eps)
  for (j in seq_len(10)) {
    step <- (digamma(a) - y) / trigamma(a)
    candidate <- a - step
    bad <- !is.finite(candidate) | candidate <= 0
    candidate[bad] <- a[bad] / 2
    a <- candidate
  }
  a
}

inv_profile_score_direct <- function(t) {
  a <- (3 + sqrt(9 + 12 * t)) / (12 * t)
  for (j in seq_len(10)) {
    step <- (profile_score_direct(a) - t) / profile_score_prime_direct(a)
    candidate <- a - step
    bad <- !is.finite(candidate) | candidate <= 0
    candidate[bad] <- a[bad] / 2
    a <- candidate
  }
  a
}

inv_profile_score <- function(t) {
  a <- (3 + sqrt(9 + 12 * t)) / (12 * t)
  for (j in seq_len(15)) {
    step <- (profile_score(a) - t) / profile_score_prime(a)
    candidate <- a - step
    bad <- !is.finite(candidate) | candidate <= 0
    candidate[bad] <- a[bad] / 2
    a <- candidate
  }
  a
}

g_dispersion <- function(phi) {
  profile_score(1 / phi)
}

g_dispersion_prime <- function(phi) {
  -profile_score_prime(1 / phi) / phi^2
}

g_dispersion_direct <- function(phi) {
  profile_score_direct(1 / phi)
}

g_dispersion_prime_direct <- function(phi) {
  -profile_score_prime_direct(1 / phi) / phi^2
}

solve_dispersion <- function(target, start) {
  phi <- pmax(start, .Machine$double.eps)
  for (j in seq_len(12)) {
    step <- (g_dispersion(phi) - target) / g_dispersion_prime(phi)
    candidate <- phi - step
    bad <- !is.finite(candidate) | candidate <= 0
    candidate[bad] <- phi[bad] / 2
    phi <- candidate
  }
  phi
}

solve_dispersion_direct <- function(target, start) {
  phi <- pmax(start, .Machine$double.eps)
  for (j in seq_len(12)) {
    step <- (g_dispersion_direct(phi) - target) / g_dispersion_prime_direct(phi)
    candidate <- phi - step
    bad <- !is.finite(candidate) | candidate <= 0
    candidate[bad] <- phi[bad] / 2
    phi <- candidate
  }
  phi
}

rmse <- function(x, truth) sqrt(mean((x - truth)^2))

# Return one row of the summaries used in the discrepancy and RMSE figures.
summarize_discrepancy <- function(values, model, parameter, method, n, growth) {
  q <- quantile(values, c(0.05, 0.50, 0.95), na.rm = TRUE)
  data.frame(
    model = model,
    parameter = parameter,
    method = method,
    n = n,
    growth = growth,
    q05 = unname(q[1]),
    median = unname(q[2]),
    q95 = unname(q[3])
  )
}

summarize_rmse <- function(estimator, exact, truth, model, parameter,
                           method, n, growth) {
  data.frame(
    model = model,
    parameter = parameter,
    method = method,
    n = n,
    growth = growth,
    rmse = rmse(estimator, truth),
    exact_rmse = rmse(exact, truth),
    ratio = rmse(estimator, truth) / rmse(exact, truth)
  )
}

# Simulate the two iid gamma designs at one (n, gamma) grid point.
simulate_iid_cell <- function(n, growth, B, seed) {
  set.seed(seed)
  alpha <- n^growth

  known_tilde <- known_hat <- numeric(B)
  known_target <- numeric(B)
  unknown_a_tilde <- unknown_a_hat <- unknown_a_hat_direct <- numeric(B)
  unknown_b_tilde <- unknown_b_hat <- numeric(B)
  profile_target <- numeric(B)
  yc_a <- yc_b <- numeric(B)

  first <- 1L
  while (first <= B) {
    b <- min(chunk_size, B - first + 1L)
    ii <- first:(first + b - 1L)
    x <- matrix(
      rgamma(b * n, shape = alpha, rate = beta_rate),
      nrow = b,
      ncol = n
    )
    log_x <- log(x)
    xbar <- rowMeans(x)
    log_g <- rowMeans(log_x)

    known_target[ii] <- log(beta_rate) + log_g
    known_tilde[ii] <- 0.5 + beta_rate * exp(log_g)
    known_hat[ii] <- inv_digamma(known_target[ii])

    profile_t <- log(xbar) - log_g
    profile_target[ii] <- profile_t
    unknown_a_tilde[ii] <- (3 + sqrt(9 + 12 * profile_t)) / (12 * profile_t)
    unknown_a_hat[ii] <- inv_profile_score(profile_t)
    unknown_a_hat_direct[ii] <- inv_profile_score_direct(profile_t)
    unknown_b_tilde[ii] <- unknown_a_tilde[ii] / xbar
    unknown_b_hat[ii] <- unknown_a_hat[ii] / xbar

    mixed_moment <- rowMeans(x * log_x) - xbar * log_g
    yc_correction <- (n - 1) / (n + 2)
    yc_a[ii] <- yc_correction * xbar / mixed_moment
    yc_b[ii] <- yc_correction / mixed_moment
    first <- first + b
  }

  discrepancy <- bind_rows(
    summarize_discrepancy(
      abs(known_tilde - known_hat) / sqrt(alpha / n),
      "Known rate", "Shape", "Approximate MLE", n, growth
    ),
    summarize_discrepancy(
      abs(unknown_a_tilde - unknown_a_hat) / (sqrt(2) * alpha / sqrt(n)),
      "Unknown rate", "Shape", "Approximate MLE", n, growth
    ),
    summarize_discrepancy(
      abs(unknown_b_tilde - unknown_b_hat) / (sqrt(2) * beta_rate / sqrt(n)),
      "Unknown rate", "Rate", "Approximate MLE", n, growth
    )
  )

  rmse_summary <- bind_rows(
    summarize_rmse(
      known_tilde, known_hat, alpha,
      "Known rate", "Shape", "Approximate MLE", n, growth
    ),
    summarize_rmse(
      unknown_a_tilde, unknown_a_hat, alpha,
      "Unknown rate", "Shape", "Approximate MLE", n, growth
    ),
    summarize_rmse(
      yc_a, unknown_a_hat, alpha,
      "Unknown rate", "Shape", "Ye-Chen (bias-corrected)", n, growth
    ),
    summarize_rmse(
      unknown_b_tilde, unknown_b_hat, beta_rate,
      "Unknown rate", "Rate", "Approximate MLE", n, growth
    ),
    summarize_rmse(
      yc_b, unknown_b_hat, beta_rate,
      "Unknown rate", "Rate", "Ye-Chen (bias-corrected)", n, growth
    )
  )

  diagnostics <- NULL
  if (n == diagnostic_n && abs(growth - diagnostic_growth) < 1e-12) {
    diagnostics <- list(
      known = data.frame(
        estimator = rep(c("Exact MLE", "Approximate MLE"), each = B),
        z = c(
          (known_hat - alpha) / sqrt(alpha / n),
          (known_tilde - alpha) / sqrt(alpha / n)
        )
      ),
      unknown = data.frame(
        z_alpha = sqrt(n) * (unknown_a_tilde - alpha) / alpha,
        z_beta = sqrt(n) * (unknown_b_tilde - beta_rate)
      )
    )
  }

  root_difference <- abs(unknown_a_hat - unknown_a_hat_direct)
  worst <- which.max(root_difference)
  numerical_accuracy <- data.frame(
    model = "Unknown rate",
    n = n,
    growth = growth,
    max_abs_stable_score_residual = max(
      abs(profile_score(unknown_a_hat) - profile_target)
    ),
    max_abs_direct_score_residual = max(
      abs(profile_score(unknown_a_hat_direct) - profile_target)
    ),
    max_abs_direct_native_residual = max(
      abs(profile_score_direct(unknown_a_hat_direct) - profile_target)
    ),
    max_abs_root_difference = max(root_difference),
    max_scaled_root_difference = max(
      root_difference / (sqrt(2) * alpha / sqrt(n))
    ),
    worst_profile_target = profile_target[worst],
    stable_root = unknown_a_hat[worst],
    direct_root = unknown_a_hat_direct[worst],
    max_abs_known_rate_residual = max(abs(digamma(known_hat) - known_target))
  )

  list(
    discrepancy = discrepancy,
    rmse = rmse_summary,
    diagnostics = diagnostics,
    numerical_accuracy = numerical_accuracy
  )
}

# Simulate the gamma-regression design and compute the raw, bias-corrected, and
# exact likelihood estimators at one (n, gamma) grid point.
simulate_regression_cell <- function(n, growth, B, seed) {
  set.seed(seed)
  phi <- n^(-growth)
  x <- seq(-1, 1, length.out = n)
  x <- x - mean(x)
  x2_sum <- sum(x^2)
  s2 <- mean(x^2)
  mu <- exp(beta_reg[1] + beta_reg[2] * x)

  beta0_raw <- beta1_raw <- phi_raw <- numeric(B)
  beta0_bc <- beta1_bc <- phi_bc <- numeric(B)
  beta0_hat <- beta1_hat <- phi_hat <- phi_hat_direct <- numeric(B)
  dispersion_target <- numeric(B)

  first <- 1L
  while (first <= B) {
    b <- min(chunk_size, B - first + 1L)
    ii <- first:(first + b - 1L)
    z <- matrix(
      rgamma(b * n, shape = 1 / phi, scale = phi),
      nrow = b,
      ncol = n
    )
    y <- sweep(z, 2, mu, `*`)
    log_y <- log(y)

    beta0_raw[ii] <- rowMeans(log_y)
    beta1_raw[ii] <- as.vector(log_y %*% x) / x2_sum
    fitted_log <- outer(beta0_raw[ii], rep(1, n)) + outer(beta1_raw[ii], x)
    phi_raw[ii] <- rowMeans((log_y - fitted_log)^2)

    v_check <- phi_raw[ii] / (1 - 2 / n)
    phi_bc[ii] <- v_check - 0.5 * v_check^2
    beta0_bc[ii] <- beta0_raw[ii] + 0.5 * v_check - v_check^2 / 6
    beta1_bc[ii] <- beta1_raw[ii]

    slope <- beta1_raw[ii]
    for (j in seq_len(10)) {
      weights <- y * exp(-outer(slope, x))
      score <- as.vector(weights %*% x) / n
      information <- as.vector(weights %*% (x^2)) / n
      slope <- slope + score / information
    }
    weights <- y * exp(-outer(slope, x))
    intercept <- log(rowMeans(weights))
    fitted_z <- weights / rowMeans(weights)
    profile_target <- -rowMeans(log(fitted_z))
    dispersion <- solve_dispersion(profile_target, pmax(phi_bc[ii], phi / 10))
    dispersion_direct <- solve_dispersion_direct(
      profile_target, pmax(phi_bc[ii], phi / 10)
    )

    beta0_hat[ii] <- intercept
    beta1_hat[ii] <- slope
    phi_hat[ii] <- dispersion
    phi_hat_direct[ii] <- dispersion_direct
    dispersion_target[ii] <- profile_target
    first <- first + b
  }

  se_beta0 <- sqrt(phi / n)
  se_beta1 <- sqrt(phi / (n * s2))
  se_phi <- sqrt(2) * phi / sqrt(n)

  discrepancy <- bind_rows(
    summarize_discrepancy(
      abs(beta0_raw - beta0_hat) / se_beta0,
      "Gamma regression", "Intercept", "Raw explicit", n, growth
    ),
    summarize_discrepancy(
      abs(beta0_bc - beta0_hat) / se_beta0,
      "Gamma regression", "Intercept", "Bias-corrected", n, growth
    ),
    summarize_discrepancy(
      abs(beta1_raw - beta1_hat) / se_beta1,
      "Gamma regression", "Slope", "Raw explicit", n, growth
    ),
    summarize_discrepancy(
      abs(phi_raw - phi_hat) / se_phi,
      "Gamma regression", "Dispersion", "Raw explicit", n, growth
    ),
    summarize_discrepancy(
      abs(phi_bc - phi_hat) / se_phi,
      "Gamma regression", "Dispersion", "Bias-corrected", n, growth
    )
  )

  rmse_summary <- bind_rows(
    summarize_rmse(
      beta0_raw, beta0_hat, beta_reg[1],
      "Gamma regression", "Intercept", "Raw explicit", n, growth
    ),
    summarize_rmse(
      beta0_bc, beta0_hat, beta_reg[1],
      "Gamma regression", "Intercept", "Bias-corrected", n, growth
    ),
    summarize_rmse(
      beta1_raw, beta1_hat, beta_reg[2],
      "Gamma regression", "Slope", "Raw explicit", n, growth
    ),
    summarize_rmse(
      phi_raw, phi_hat, phi,
      "Gamma regression", "Dispersion", "Raw explicit", n, growth
    ),
    summarize_rmse(
      phi_bc, phi_hat, phi,
      "Gamma regression", "Dispersion", "Bias-corrected", n, growth
    )
  )

  diagnostics <- NULL
  if (n == diagnostic_n && abs(growth - diagnostic_growth) < 1e-12) {
    diagnostics <- bind_rows(
      data.frame(
        parameter = "Intercept",
        exact = sqrt(n / phi) * (beta0_hat - beta_reg[1]),
        corrected = sqrt(n / phi) * (beta0_bc - beta_reg[1])
      ),
      data.frame(
        parameter = "Slope",
        exact = sqrt(n * s2 / phi) * (beta1_hat - beta_reg[2]),
        corrected = sqrt(n * s2 / phi) * (beta1_bc - beta_reg[2])
      ),
      data.frame(
        parameter = "Dispersion",
        exact = sqrt(n / 2) * (phi_hat - phi) / phi,
        corrected = sqrt(n / 2) * (phi_bc - phi) / phi
      )
    )
  }

  root_difference <- abs(phi_hat - phi_hat_direct)
  worst <- which.max(root_difference)
  numerical_accuracy <- data.frame(
    model = "Gamma regression",
    n = n,
    growth = growth,
    max_abs_stable_score_residual = max(
      abs(g_dispersion(phi_hat) - dispersion_target)
    ),
    max_abs_direct_score_residual = max(
      abs(g_dispersion(phi_hat_direct) - dispersion_target)
    ),
    max_abs_direct_native_residual = max(
      abs(g_dispersion_direct(phi_hat_direct) - dispersion_target)
    ),
    max_abs_root_difference = max(root_difference),
    max_scaled_root_difference = max(
      root_difference / (sqrt(2) * phi / sqrt(n))
    ),
    worst_profile_target = dispersion_target[worst],
    stable_root = phi_hat[worst],
    direct_root = phi_hat_direct[worst],
    max_abs_known_rate_residual = NA_real_
  )

  list(
    discrepancy = discrepancy,
    rmse = rmse_summary,
    diagnostics = diagnostics,
    numerical_accuracy = numerical_accuracy
  )
}

# Each grid point has a deterministic seed, so results are reproducible even if
# cells are run separately or in a different order.
message(
  "Unified gamma study: ", length(n_grid) * length(growth_grid),
  " cells per sampling design; B = ", B
)

iid_results <- vector("list", length(n_grid) * length(growth_grid))
reg_results <- vector("list", length(n_grid) * length(growth_grid))
cell <- 1L
for (g_index in seq_along(growth_grid)) {
  growth <- growth_grid[g_index]
  for (n_index in seq_along(n_grid)) {
    n <- n_grid[n_index]
    message("  gamma = ", growth, ", n = ", n)
    iid_results[[cell]] <- simulate_iid_cell(
      n, growth, B,
      seed = 310000L + 1000L * g_index + n_index
    )
    reg_results[[cell]] <- simulate_regression_cell(
      n, growth, B,
      seed = 710000L + 1000L * g_index + n_index
    )
    cell <- cell + 1L
  }
}

discrepancy_summary <- bind_rows(
  lapply(iid_results, `[[`, "discrepancy"),
  lapply(reg_results, `[[`, "discrepancy")
)
rmse_summary <- bind_rows(
  lapply(iid_results, `[[`, "rmse"),
  lapply(reg_results, `[[`, "rmse")
)
numerical_accuracy <- bind_rows(
  lapply(iid_results, `[[`, "numerical_accuracy"),
  lapply(reg_results, `[[`, "numerical_accuracy")
)

iid_diag <- iid_results[[
  which(vapply(iid_results, function(x) !is.null(x$diagnostics), logical(1)))
]]$diagnostics
reg_diag <- reg_results[[
  which(vapply(reg_results, function(x) !is.null(x$diagnostics), logical(1)))
]]$diagnostics

panel_levels <- c(
  "Gamma known rate: shape",
  "Gamma unknown rate: shape",
  "Gamma unknown rate: rate",
  "Gamma regression: intercept",
  "Gamma regression: slope",
  "Gamma regression: dispersion"
)

make_panel <- function(model, parameter) {
  if (model == "Known rate") return("Gamma known rate: shape")
  if (model == "Unknown rate" && parameter == "Shape") return("Gamma unknown rate: shape")
  if (model == "Unknown rate" && parameter == "Rate") return("Gamma unknown rate: rate")
  paste0("Gamma regression: ", tolower(parameter))
}

discrepancy_summary$panel <- factor(
  mapply(make_panel, discrepancy_summary$model, discrepancy_summary$parameter),
  levels = panel_levels
)
rmse_summary$panel <- factor(
  mapply(make_panel, rmse_summary$model, rmse_summary$parameter),
  levels = panel_levels
)
discrepancy_summary$growth_factor <- factor(discrepancy_summary$growth)
rmse_summary$growth_factor <- factor(rmse_summary$growth)

method_levels <- c(
  "Approximate MLE", "Ye-Chen (bias-corrected)", "Raw explicit", "Bias-corrected"
)
method_linetypes <- c(
  "Approximate MLE" = "solid",
  "Ye-Chen (bias-corrected)" = "longdash",
  "Raw explicit" = "dotted",
  "Bias-corrected" = "solid"
)
method_shapes <- c(
  "Approximate MLE" = 16,
  "Ye-Chen (bias-corrected)" = 17,
  "Raw explicit" = 1,
  "Bias-corrected" = 15
)

p_discrepancy <- ggplot(
  discrepancy_summary,
  aes(
    x = n,
    y = median,
    colour = growth_factor,
    linetype = method,
    shape = method,
    group = interaction(growth_factor, method)
  )
) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.35, linetype = "longdash") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.25, stroke = 0.45) +
  facet_wrap(vars(panel), ncol = 3, scales = "free_y") +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(
    breaks = 10^seq(-12, 4, by = 1),
    labels = power_of_ten_labels
  ) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  scale_linetype_manual(values = method_linetypes, breaks = method_levels, name = NULL) +
  scale_shape_manual(values = method_shapes, breaks = method_levels, name = NULL) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)
  ) +
  labs(
    x = expression(italic(n)),
    y = "Median approximate-exact discrepancy\n(asymptotic standard error units)"
  ) +
  paper_theme

p_rmse <- ggplot(
  rmse_summary,
  aes(
    x = n,
    y = ratio,
    colour = growth_factor,
    linetype = method,
    shape = method,
    group = interaction(growth_factor, method)
  )
) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.4, linetype = "longdash") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.25, stroke = 0.45) +
  facet_wrap(vars(panel), ncol = 3, scales = "free_y") +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  scale_linetype_manual(values = method_linetypes, breaks = method_levels, name = NULL) +
  scale_shape_manual(values = method_shapes, breaks = method_levels, name = NULL) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)
  ) +
  labs(
    x = expression(italic(n)),
    y = "RMSE relative to the exact MLE"
  ) +
  paper_theme

# Cairo creates figures with the manuscript font. Ghostscript rewrites each
# temporary file with embedded fonts for reliable inclusion in LaTeX.
cm_cairo_pdf <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(
    filename, width = width, height = height,
    family = figure_family, symbolfamily = figure_family, ...
  )
}

convert_pdf_for_latex <- function(input_pdf, output_pdf) {
  gs <- Sys.which("gs")
  if (!nzchar(gs)) stop("Ghostscript is required to create the final figure PDFs.")
  output <- system2(
    gs,
    c(
      "-dNOPAUSE", "-dBATCH", "-dSAFER",
      "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.5",
      "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true",
      "-dSubsetFonts=true",
      paste0("-sOutputFile=", shQuote(output_pdf)),
      shQuote(input_pdf)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop(paste(output, collapse = "\n"))
  }
}

save_ggplot <- function(plot, stem, width, height) {
  final_pdf <- file.path(fig_dir, paste0(stem, ".pdf"))
  raw_pdf <- tempfile(fileext = ".pdf")
  ggsave(
    raw_pdf, plot,
    width = width, height = height, device = cm_cairo_pdf
  )
  convert_pdf_for_latex(raw_pdf, final_pdf)
  unlink(raw_pdf)
}

save_ggplot(p_discrepancy, "fig01_unified_discrepancy", 8.1, 5.25)
save_ggplot(p_rmse, "fig02_unified_rmse", 8.1, 5.25)

p_known_diag <- ggplot(iid_diag$known, aes(x = z, colour = estimator, linetype = estimator)) +
  geom_density(linewidth = 0.75, adjust = 1.1, key_glyph = "path") +
  stat_function(
    fun = dnorm,
    aes(colour = "Standard normal", linetype = "Standard normal"),
    linewidth = 0.65,
    key_glyph = "path"
  ) +
  scale_colour_manual(
    values = c("Exact MLE" = blue, "Approximate MLE" = orange, "Standard normal" = grey),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Exact MLE" = "solid", "Approximate MLE" = "dashed", "Standard normal" = "dotted"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-4, 4)) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x = "Standardized shape error",
    y = "Density",
    title = "(a) Gamma known rate: marginal limit"
  ) +
  paper_theme +
  theme(legend.position = "bottom", plot.title = element_text(size = 11, face = "plain"))

p_unknown_diag <- ggplot(iid_diag$unknown, aes(x = z_alpha, y = z_beta)) +
  geom_point(colour = blue, alpha = 0.18, size = 0.6) +
  geom_abline(intercept = 0, slope = beta_rate, linetype = "dashed", colour = grey) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-10, 10)) +
  labs(
    x = expression(sqrt(italic(n)) * (tilde(alpha)[n] - alpha[n]) / alpha[n]),
    y = expression(sqrt(italic(n)) * (tilde(beta)[n] - beta)),
    title = "(b) Gamma unknown rate: joint limit"
  ) +
  paper_theme +
  theme(legend.position = "none", plot.title = element_text(size = 11, face = "plain"))

qq_probs <- ppoints(B)
reg_qq <- bind_rows(lapply(split(reg_diag, reg_diag$parameter), function(d) {
  data.frame(
    parameter = d$parameter[1],
    theoretical = qnorm(qq_probs),
    `Exact MLE` = sort(d$exact),
    `Bias-corrected` = sort(d$corrected),
    check.names = FALSE
  )
})) %>%
  tidyr::pivot_longer(
    cols = c("Exact MLE", "Bias-corrected"),
    names_to = "estimator",
    values_to = "empirical"
  )
reg_qq$parameter <- factor(reg_qq$parameter, levels = c("Intercept", "Slope", "Dispersion"))

p_reg_diag <- ggplot(reg_qq, aes(theoretical, empirical, colour = estimator, linetype = estimator)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.65) +
  facet_wrap(vars(parameter), nrow = 1) +
  scale_colour_manual(values = c("Exact MLE" = blue, "Bias-corrected" = green), name = NULL) +
  scale_linetype_manual(values = c("Exact MLE" = "solid", "Bias-corrected" = "dashed"), name = NULL) +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5)) +
  labs(
    x = "Standard normal quantile",
    y = "Empirical quantile",
    title = "(c) Gamma regression: standardized marginal limits"
  ) +
  paper_theme +
  theme(plot.title = element_text(size = 11, face = "plain"))

# Arrange the known-rate, unknown-rate, and regression diagnostics in the same
# three-part layout used in the manuscript.
draw_diagnostics <- function() {
  grid::grid.newpage()
  known_legend_table <- ggplotGrob(p_known_diag)
  legend_index <- which(vapply(
    known_legend_table$grobs,
    function(grob) identical(grob$name, "guide-box"),
    logical(1)
  ))
  known_legend <- known_legend_table$grobs[[legend_index]]

  known_plot <- ggplotGrob(p_known_diag + theme(legend.position = "none"))
  unknown_plot <- ggplotGrob(p_unknown_diag)
  aligned_widths <- grid::unit.pmax(known_plot$widths, unknown_plot$widths)
  aligned_heights <- grid::unit.pmax(known_plot$heights, unknown_plot$heights)
  known_plot$widths <- aligned_widths
  unknown_plot$widths <- aligned_widths
  known_plot$heights <- aligned_heights
  unknown_plot$heights <- aligned_heights

  layout <- grid::grid.layout(
    nrow = 3,
    ncol = 2,
    heights = grid::unit(c(0.43, 0.08, 0.49), "null"),
    widths = grid::unit(c(0.5, 0.5), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(known_plot)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.draw(unknown_plot)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.draw(known_legend)
  grid::upViewport()
  print(p_reg_diag, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2))
  grid::upViewport()
}

diagnostic_raw_pdf <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(
  diagnostic_raw_pdf,
  width = 8.1, height = 6.6,
  family = figure_family, symbolfamily = figure_family
)
draw_diagnostics()
grDevices::dev.off()
convert_pdf_for_latex(
  diagnostic_raw_pdf,
  file.path(fig_dir, "fig03_unified_diagnostics.pdf")
)
unlink(diagnostic_raw_pdf)

# Save only the compact tables needed to reproduce the figures and the
# high-precision numerical audit.
write.csv(discrepancy_summary, file.path(out_dir, "discrepancy_summary.csv"), row.names = FALSE)
write.csv(rmse_summary, file.path(out_dir, "rmse_summary.csv"), row.names = FALSE)
write.csv(
  numerical_accuracy,
  file.path(out_dir, "numerical_accuracy.csv"),
  row.names = FALSE
)

config <- c(
  "iid_seed=310000 + 1000*growth_index + n_index",
  "regression_seed=710000 + 1000*growth_index + n_index",
  paste0("B=", B),
  paste0("chunk_size=", chunk_size),
  paste0("n_grid=", paste(n_grid, collapse = ",")),
  paste0("growth_grid=", paste(growth_grid, collapse = ",")),
  paste0("iid_rate=", beta_rate),
  paste0("regression_beta=", paste(beta_reg, collapse = ",")),
  "regression_design=(1,x), with x equally spaced and centered on [-1,1]",
  paste0("diagnostic_cell=(n,gamma)=", diagnostic_n, ",", diagnostic_growth),
  "profile_score=direct below shape 50; cancellation-free Bernoulli expansion at and above 50",
  "numerical_audit=stable roots compared with the direct special-function formulation"
)
writeLines(config, file.path(out_dir, "simulation_config.txt"))
message("Done. Outputs written to ", out_dir, " and ", fig_dir)
