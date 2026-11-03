#!/usr/bin/env Rscript

# Monte Carlo experiments for the generalized Pareto example in the paper.
#
# The script computes the restricted global MLE, the constrained surrogate-root
# estimator, and the Hosking--Wallis probability-weighted-moment benchmark. It
# produces the three GPD figures and the summaries quoted in the manuscript.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

options(scipen = 4)

# Output locations can be overridden without changing the script.
out_dir <- Sys.getenv("GPD_SIM_OUT_DIR", "output/summaries/gpd")
fig_dir <- Sys.getenv("GPD_SIM_FIG_DIR", "output/figures/gpd")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("GPD_SIM_B", "2000"))
chunk_size <- as.integer(Sys.getenv("GPD_SIM_CHUNK", "100"))
quick <- identical(Sys.getenv("GPD_SIM_QUICK", "0"), "1")

if (quick) {
  n_grid <- c(250L, 1000L)
} else {
  n_grid <- c(250L, 500L, 1000L, 2500L, 5000L, 10000L, 25000L)
}
growth_grid <- c(0.20, 0.30, 0.55)
shape_coefficient <- as.numeric(Sys.getenv("GPD_SIM_SHAPE_COEFFICIENT", "0.30"))
shape_upper_bound <- 0.16
diagnostic_n <- max(n_grid)
diagnostic_growth <- 0.30
figure_family <- Sys.getenv("FIGURE_FONT", "CMU Serif")

blue <- "#0072B2"
orange <- "#E69F00"
purple <- "#CC79A7"
grey <- "#4D4D4D"
growth_colours <- c(
  "0.2" = blue,
  "0.3" = orange,
  "0.55" = purple
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

rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))

# Solve the quadratic surrogate equations. If the raw root is not in the
# restricted parameter space, use the shifted-exponential boundary estimator.
gpd_approximate_estimator <- function(m1, m2, m3, shape_upper_bound) {
  Cn <- 6 * m1 * m2 + 8 * m3
  En <- 8 * m1 * m3 - 3 * m2^2
  discriminant <- Cn^2 - 48 * m2 * En
  sigma_raw <- xi_raw <- NA_real_
  if (!is.finite(discriminant) || discriminant < 0 ||
      !is.finite(m2) || m2 <= 0) {
    return(c(
      sigma = m1, xi = 0, interior = 0,
      sigma_raw = sigma_raw, xi_raw = xi_raw
    ))
  }
  sigma_raw <- (Cn - sqrt(discriminant)) / (24 * m2)
  xi_raw <- 2 * sigma_raw * (m1 - sigma_raw) / m2
  interior <- is.finite(sigma_raw) && is.finite(xi_raw) && sigma_raw > 0 &&
    xi_raw >= 0 && xi_raw <= shape_upper_bound
  sigma <- if (interior) sigma_raw else m1
  xi <- if (interior) xi_raw else 0
  c(
    sigma = sigma, xi = xi, interior = as.integer(interior),
    sigma_raw = sigma_raw, xi_raw = xi_raw
  )
}

# Hosking--Wallis probability-weighted-moment estimator applied to the excesses
# above the sample minimum.
gpd_pwm_estimator <- function(y, shape_upper_bound) {
  n <- length(y)
  ordered <- sort(y)
  a0 <- mean(ordered)
  a1 <- sum(ordered * (n - seq_len(n)) / (n - 1)) / n
  denominator <- a0 - 2 * a1
  xi <- 2 - a0 / denominator
  sigma <- 2 * a0 * a1 / denominator
  valid <- is.finite(sigma) && is.finite(xi) && sigma > 0 &&
    xi >= 0 && xi <= shape_upper_bound
  c(sigma = sigma, xi = xi, valid = as.integer(valid))
}

# Restricted two-parameter log-likelihood for excesses above the location
# estimate. The xi=0 branch is evaluated by its exponential limit.
gpd_loglik <- function(y, sigma, xi) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(xi) || xi < 0) return(-Inf)
  n <- length(y)
  if (xi == 0) return(-n * log(sigma) - sum(y) / sigma)
  z <- xi * y / sigma
  if (any(z <= -1) || any(!is.finite(z))) return(-Inf)
  -n * log(sigma) - (1 + 1 / xi) * sum(log1p(z))
}

# Root-preserving likelihood equations and their Jacobian for xi>0.
gpd_system <- function(y, sigma, xi) {
  t <- y / sigma
  d <- 1 + xi * t
  Bn <- mean(log1p(xi * t) / xi)
  An <- mean((1 + xi) * t / d)
  g1 <- Bn - 1
  g2 <- (Bn - An) / xi
  Bs <- -mean(t / (sigma * d))
  Be <- mean((xi * t / d - log1p(xi * t)) / xi^2)
  As <- -mean((1 + xi) * t / (sigma * d^2))
  Ae <- mean(t / d - (1 + xi) * t^2 / d^2)
  J <- matrix(
    c(Bs, Be, (Bs - As) / xi, (Be - Ae - g2) / xi),
    nrow = 2,
    byrow = TRUE
  )
  list(g = c(g1, g2), J = J, A = An)
}

gpd_profile_scale <- function(y, xi, start) {
  if (xi == 0) return(mean(y))
  sigma <- max(start, mean(y) * 1e-6, .Machine$double.eps)
  score <- function(value) {
    t <- y / value
    mean((1 + xi) * t / (1 + xi * t)) - 1
  }
  for (iteration in seq_len(25L)) {
    t <- y / sigma
    d <- 1 + xi * t
    value <- mean((1 + xi) * t / d) - 1
    derivative <- -mean((1 + xi) * t / (sigma * d^2))
    if (abs(value) < 1e-12) return(sigma)
    candidate <- sigma - value / derivative
    if (!is.finite(candidate) || candidate <= 0) candidate <- sigma / 2
    if (abs(score(candidate)) > abs(value)) candidate <- (candidate + sigma) / 2
    sigma <- candidate
  }
  lower <- max(mean(y) * 1e-12, .Machine$double.eps)
  upper <- max(mean(y), sigma, 1)
  while (score(upper) > 0 && upper < 1e12) upper <- 2 * upper
  for (iteration in seq_len(70L)) {
    middle <- (lower + upper) / 2
    if (score(middle) > 0) lower <- middle else upper <- middle
  }
  (lower + upper) / 2
}

# Damped Newton solver for an interior likelihood root.
gpd_interior_newton <- function(y, start, shape_upper_bound) {
  sigma <- max(start[1], mean(y) * 1e-5, .Machine$double.eps)
  xi <- min(max(start[2], 1e-7), shape_upper_bound - 1e-7)
  current <- gpd_system(y, sigma, xi)
  residual <- max(abs(current$g))
  for (iteration in seq_len(35L)) {
    if (!all(is.finite(current$g)) || !all(is.finite(current$J))) break
    if (residual < 1e-9) {
      return(list(
        sigma = sigma, xi = xi, residual = residual,
        converged = TRUE, iterations = iteration
      ))
    }
    step <- tryCatch(
      solve(current$J, current$g),
      error = function(e) c(NA_real_, NA_real_)
    )
    if (any(!is.finite(step))) break
    step_length <- 1
    accepted <- FALSE
    for (halving in seq_len(24L)) {
      candidate_sigma <- sigma - step_length * step[1]
      candidate_xi <- xi - step_length * step[2]
      if (is.finite(candidate_sigma) && candidate_sigma > 0 &&
          is.finite(candidate_xi) && candidate_xi > 1e-8 &&
          candidate_xi < shape_upper_bound - 1e-8) {
        candidate <- gpd_system(y, candidate_sigma, candidate_xi)
        candidate_residual <- max(abs(candidate$g))
        if (is.finite(candidate_residual) && candidate_residual < residual) {
          sigma <- candidate_sigma
          xi <- candidate_xi
          current <- candidate
          residual <- candidate_residual
          accepted <- TRUE
          break
        }
      }
      step_length <- step_length / 2
    }
    if (!accepted) break
  }
  list(
    sigma = sigma, xi = xi, residual = residual,
    converged = residual < 1e-7, iterations = 35L
  )
}

# Compare all converged interior roots with the xi=0 and xi=shape_upper_bound
# boundary solutions, then return the candidate with the largest likelihood.
gpd_exact_mle <- function(y, shape_upper_bound, starts) {
  m1 <- mean(y)
  m2 <- mean(y^2)
  candidates <- list()
  candidate_index <- 1L
  for (start in starts) {
    if (all(is.finite(start))) {
      answer <- gpd_interior_newton(y, start, shape_upper_bound)
      if (isTRUE(answer$converged)) {
        candidates[[candidate_index]] <- c(
          sigma = answer$sigma,
          xi = answer$xi,
          loglik = gpd_loglik(y, answer$sigma, answer$xi),
          residual = answer$residual,
          boundary = 0
        )
        candidate_index <- candidate_index + 1L
      }
    }
  }

  sigma_zero <- m1
  q_zero <- -m1 / sigma_zero + m2 / (2 * sigma_zero^2)
  candidates[[candidate_index]] <- c(
    sigma = sigma_zero,
    xi = 0,
    loglik = gpd_loglik(y, sigma_zero, 0),
    residual = 0,
    boundary = 1
  )
  candidate_index <- candidate_index + 1L

  sigma_upper <- gpd_profile_scale(y, shape_upper_bound, m1)
  score_upper <- gpd_system(y, sigma_upper, shape_upper_bound)
  candidates[[candidate_index]] <- c(
    sigma = sigma_upper,
    xi = shape_upper_bound,
    loglik = gpd_loglik(y, sigma_upper, shape_upper_bound),
    residual = abs(score_upper$A - 1),
    boundary = 2
  )

  candidate_matrix <- do.call(rbind, candidates)
  best <- candidate_matrix[which.max(candidate_matrix[, "loglik"]), ]
  boundary_ok <- (best["boundary"] == 0) ||
    (best["boundary"] == 1 && q_zero <= 1e-6) ||
    (best["boundary"] == 2 && score_upper$g[2] >= -1e-6)
  residual_ok <- is.finite(best["residual"]) && best["residual"] <= 1e-6
  failure <- as.integer(!is.finite(best["loglik"]) || !boundary_ok || !residual_ok)
  c(
    sigma = unname(best["sigma"]),
    xi = unname(best["xi"]),
    loglik = unname(best["loglik"]),
    residual = unname(best["residual"]),
    boundary = unname(best["boundary"]),
    failure = failure
  )
}

# Inverse-transform sampler for the GPD parameterization used in the paper.
r_gpd <- function(n, sigma, xi, mu = 0) {
  exponential <- rexp(n)
  if (xi == 0) return(mu + sigma * exponential)
  mu + sigma * expm1(xi * exponential) / xi
}

# Simulate and estimate one (n, gamma) cell in memory-bounded chunks.
simulate_cell <- function(
    n, growth, B, chunk_size, shape_coefficient, shape_upper_bound, seed) {
  set.seed(seed)
  xi_true <- shape_coefficient * n^(-growth)
  pieces <- vector("list", ceiling(B / chunk_size))
  piece_index <- 1L

  for (first in seq.int(1L, B, by = chunk_size)) {
    batch <- min(chunk_size, B - first + 1L)
    x <- matrix(r_gpd(batch * n, sigma = 1, xi = xi_true), nrow = batch)
    mu_hat <- apply(x, 1, min)
    y <- x - mu_hat
    m1 <- rowMeans(y)
    m2 <- rowMeans(y^2)
    m3 <- rowMeans(y^3)
    rows <- vector("list", batch)

    for (index in seq_len(batch)) {
      approximate <- gpd_approximate_estimator(
        m1[index], m2[index], m3[index], shape_upper_bound
      )
      pwm <- gpd_pwm_estimator(y[index, ], shape_upper_bound)
      starts <- list(
        c(approximate["sigma"], approximate["xi"]),
        c(pwm["sigma"], pwm["xi"]),
        c(m1[index], min(0.05, shape_upper_bound / 2))
      )
      exact <- gpd_exact_mle(y[index, ], shape_upper_bound, starts)
      rows[[index]] <- data.frame(
        growth = growth,
        n = n,
        xi_true = xi_true,
        mu_hat = mu_hat[index],
        sigma_exact = exact["sigma"],
        xi_exact = exact["xi"],
        exact_boundary = exact["boundary"],
        exact_residual = exact["residual"],
        exact_failure = exact["failure"],
        sigma_approximate = approximate["sigma"],
        xi_approximate = approximate["xi"],
        approximate_interior = approximate["interior"],
        sigma_approximate_raw = approximate["sigma_raw"],
        xi_approximate_raw = approximate["xi_raw"],
        sigma_pwm = pwm["sigma"],
        xi_pwm = pwm["xi"],
        pwm_valid = pwm["valid"]
      )
    }
    pieces[[piece_index]] <- bind_rows(rows)
    piece_index <- piece_index + 1L
  }

  answer <- bind_rows(pieces)
  answer$approximate_comparison_valid <-
    answer$exact_failure == 0
  answer$pwm_comparison_valid <-
    answer$exact_failure == 0 & answer$pwm_valid == 1
  answer$z_sigma_exact <- sqrt(n) * (answer$sigma_exact - 1)
  answer$z_xi_exact <- sqrt(n) * (answer$xi_exact - xi_true)
  answer$z_sigma_approximate <- sqrt(n) * (answer$sigma_approximate - 1)
  answer$z_xi_approximate <- sqrt(n) * (answer$xi_approximate - xi_true)
  answer$location_scaled <- n * answer$mu_hat
  answer
}

summarize_cell <- function(dat) {
  exact_valid <- dat$exact_failure == 0
  approximate_valid <- exact_valid
  pwm_valid <- dat$pwm_comparison_valid
  xi_true <- dat$xi_true[1]

  summarize_parameter <- function(parameter, standard_error) {
    exact <- dat[[paste0(parameter, "_exact")]]
    approximate <- dat[[paste0(parameter, "_approximate")]]
    pwm <- dat[[paste0(parameter, "_pwm")]]
    truth <- if (parameter == "sigma") 1 else xi_true
    discrepancy <- abs(approximate[approximate_valid] - exact[approximate_valid]) /
      standard_error
    data.frame(
      parameter = if (parameter == "sigma") "Scale" else "Shape",
      discrepancy_median = median(discrepancy, na.rm = TRUE),
      discrepancy_q05 = unname(quantile(discrepancy, 0.05, na.rm = TRUE)),
      discrepancy_q95 = unname(quantile(discrepancy, 0.95, na.rm = TRUE)),
      approximate_rmse_ratio =
        rmse(approximate[approximate_valid] - truth) /
        rmse(exact[approximate_valid] - truth),
      pwm_rmse_ratio =
        rmse(pwm[pwm_valid] - truth) /
        rmse(exact[pwm_valid] - truth),
      exact_bias = mean(exact[exact_valid] - truth, na.rm = TRUE),
      exact_rmse = rmse(exact[exact_valid] - truth),
      approximate_bias = mean(approximate[approximate_valid] - truth, na.rm = TRUE),
      approximate_rmse = rmse(approximate[approximate_valid] - truth),
      pwm_bias = mean(pwm[pwm_valid] - truth, na.rm = TRUE),
      pwm_rmse = rmse(pwm[pwm_valid] - truth)
    )
  }

  parameter_summary <- bind_rows(
    summarize_parameter("sigma", sqrt(2 / dat$n[1])),
    summarize_parameter("xi", 1 / sqrt(dat$n[1]))
  )
  parameter_summary$growth <- dat$growth[1]
  parameter_summary$n <- dat$n[1]
  parameter_summary$xi_true <- xi_true
  parameter_summary$replications <- nrow(dat)
  parameter_summary$exact_failure_rate <- mean(dat$exact_failure != 0)
  parameter_summary$exact_lower_boundary_rate <- mean(dat$exact_boundary == 1)
  parameter_summary$exact_upper_boundary_rate <- mean(dat$exact_boundary == 2)
  parameter_summary$raw_root_replacement_rate <- mean(dat$approximate_interior != 1)
  parameter_summary$pwm_invalid_rate <- mean(dat$pwm_valid != 1)
  parameter_summary$approximate_comparison_n <- sum(approximate_valid)
  parameter_summary$pwm_comparison_n <- sum(pwm_valid)
  parameter_summary$max_exact_residual <- max(dat$exact_residual, na.rm = TRUE)
  parameter_summary
}

# Deterministic cell-specific seeds make the results reproducible regardless of
# the order in which the grid is evaluated.
message(
  "GPD study: ", length(n_grid) * length(growth_grid),
  " cells; B = ", B
)

cells <- vector("list", length(n_grid) * length(growth_grid))
cell_index <- 1L
for (growth_index in seq_along(growth_grid)) {
  growth <- growth_grid[growth_index]
  for (n_index in seq_along(n_grid)) {
    n <- n_grid[n_index]
    message("  gamma = ", growth, ", n = ", n)
    cells[[cell_index]] <- simulate_cell(
      n, growth, B, chunk_size, shape_coefficient, shape_upper_bound,
      seed = 850000L + 1000L * growth_index + n_index
    )
    cell_index <- cell_index + 1L
  }
}

replications <- bind_rows(cells)
summary <- bind_rows(lapply(cells, summarize_cell))
summary <- summary[order(summary$growth, summary$n, summary$parameter), ]
summary$growth_factor <- factor(summary$growth)
summary$parameter <- factor(summary$parameter, levels = c("Scale", "Shape"))

discrepancy_summary <- summary %>%
  select(
    growth, growth_factor, n, parameter,
    discrepancy_median, discrepancy_q05, discrepancy_q95
  )

rmse_summary <- bind_rows(
  summary %>% transmute(
    growth, growth_factor, n, parameter,
    method = "Approximate MLE", ratio = approximate_rmse_ratio
  ),
  summary %>% transmute(
    growth, growth_factor, n, parameter,
    method = "Hosking-Wallis PWM", ratio = pwm_rmse_ratio
  )
)
rmse_summary$method <- factor(
  rmse_summary$method,
  levels = c("Approximate MLE", "Hosking-Wallis PWM")
)

event_summary <- summary %>%
  filter(parameter == "Scale") %>%
  select(
    growth, growth_factor, n,
    exact_lower_boundary_rate,
    exact_upper_boundary_rate,
    raw_root_replacement_rate
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      exact_lower_boundary_rate,
      exact_upper_boundary_rate,
      raw_root_replacement_rate
    ),
    names_to = "event",
    values_to = "rate"
  )
event_summary$event <- factor(
  event_summary$event,
  levels = c(
    "exact_lower_boundary_rate",
    "exact_upper_boundary_rate",
    "raw_root_replacement_rate"
  ),
  labels = c(
    "'Exact MLE: ' * hat(xi)[n] == 0",
    "'Exact MLE: ' * hat(xi)[n] == bar(xi)",
    "'Raw surrogate root outside ' * Theta"
  )
)

p_discrepancy <- ggplot(
  discrepancy_summary,
  aes(
    x = n, y = discrepancy_median,
    colour = growth_factor, group = growth_factor
  )
) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.35, linetype = "longdash") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.25) +
  facet_grid(. ~ parameter) +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(
    breaks = 10^seq(-3, 2, by = 1),
    labels = power_of_ten_labels
  ) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(
    x = expression(italic(n)),
    y = "Median approximate-exact discrepancy\n(asymptotic standard error units)"
  ) +
  paper_theme

method_linetypes <- c(
  "Approximate MLE" = "solid",
  "Hosking-Wallis PWM" = "longdash"
)
method_shapes <- c(
  "Approximate MLE" = 16,
  "Hosking-Wallis PWM" = 17
)

p_rmse <- ggplot(
  rmse_summary,
  aes(
    x = n, y = ratio,
    colour = growth_factor,
    linetype = method,
    shape = method,
    group = interaction(growth_factor, method)
  )
) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.4, linetype = "longdash") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.25, stroke = 0.45) +
  facet_grid(. ~ parameter) +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  scale_linetype_manual(values = method_linetypes, name = NULL) +
  scale_shape_manual(values = method_shapes, name = NULL) +
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

diagnostic <- replications[
  replications$n == diagnostic_n &
    abs(replications$growth - diagnostic_growth) < 1e-12,
]
diagnostic <- diagnostic[diagnostic$exact_failure == 0, ]
diagnostic_approximate <- diagnostic

joint_data <- bind_rows(
  data.frame(
    estimator = "Exact MLE",
    scale_error = diagnostic$z_sigma_exact,
    shape_error = diagnostic$z_xi_exact
  ),
  data.frame(
    estimator = "Approximate MLE",
    scale_error = diagnostic_approximate$z_sigma_approximate,
    shape_error = diagnostic_approximate$z_xi_approximate
  )
)
joint_data$estimator <- factor(
  joint_data$estimator,
  levels = c("Exact MLE", "Approximate MLE")
)

ellipse_points <- function(Sigma, level = 0.90, points = 300) {
  decomposition <- eigen(Sigma, symmetric = TRUE)
  angle <- seq(0, 2 * pi, length.out = points)
  radius <- sqrt(qchisq(level, df = 2))
  coordinates <- radius * decomposition$vectors %*%
    diag(sqrt(decomposition$values)) %*%
    rbind(cos(angle), sin(angle))
  data.frame(scale_error = coordinates[1, ], shape_error = coordinates[2, ])
}

target_ellipse <- ellipse_points(matrix(c(2, -1, -1, 1), nrow = 2))

p_joint <- ggplot(
  joint_data,
  aes(x = scale_error, y = shape_error, colour = estimator, shape = estimator)
) +
  geom_point(alpha = 0.16, size = 0.65) +
  geom_path(
    data = target_ellipse,
    aes(x = scale_error, y = shape_error),
    inherit.aes = FALSE,
    colour = grey,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_colour_manual(
    values = c("Exact MLE" = blue, "Approximate MLE" = orange),
    name = NULL
  ) +
  scale_shape_manual(values = c("Exact MLE" = 16, "Approximate MLE" = 17), name = NULL) +
  coord_equal(xlim = c(-5, 5), ylim = c(-5, 5)) +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1, size = 1.5), nrow = 1),
    shape = "none"
  ) +
  labs(
    x = "Scaled scale error",
    y = "Scaled shape error",
    title = "(a) GPD: joint limit"
  ) +
  paper_theme +
  theme(plot.title = element_text(size = 11, face = "plain"))

location_quantiles <- data.frame(
  theoretical = qexp(ppoints(nrow(diagnostic))),
  empirical = sort(diagnostic$location_scaled)
)
p_location <- ggplot(location_quantiles, aes(x = theoretical, y = empirical)) +
  geom_abline(intercept = 0, slope = 1, colour = grey, linetype = "dashed", linewidth = 0.6) +
  geom_point(colour = blue, alpha = 0.28, size = 0.75) +
  coord_cartesian(xlim = c(0, 7.5), ylim = c(0, 7.5)) +
  labs(
    x = expression(Exp(1)~"quantiles"),
    y = expression("Quantiles of " * italic(n) * (hat(mu)[n] - mu) / sigma),
    title = "(b) GPD: location limit"
  ) +
  paper_theme +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, face = "plain")
  )

p_events <- ggplot(
  event_summary,
  aes(x = n, y = rate, colour = growth_factor, group = growth_factor)
) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.15) +
  facet_grid(. ~ event, labeller = label_parsed) +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(
    x = expression(italic(n)),
    y = "Empirical frequency",
    title = "(c) GPD: boundary and localization events"
  ) +
  paper_theme +
  theme(plot.title = element_text(size = 11, face = "plain"))

# Cairo and Ghostscript produce PDF figures with embedded manuscript fonts.
cm_cairo_pdf <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(
    filename, width = width, height = height,
    family = figure_family, symbolfamily = figure_family, ...
  )
}

convert_pdf_for_latex <- function(input_pdf, output_pdf) {
  ghostscript <- Sys.which("gs")
  if (!nzchar(ghostscript)) stop("Ghostscript is required to create the final figure PDFs.")
  output <- system2(
    ghostscript,
    c(
      "-dNOPAUSE", "-dBATCH", "-dSAFER",
      "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.5",
      "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true", "-dSubsetFonts=true",
      paste0("-sOutputFile=", shQuote(output_pdf)),
      shQuote(input_pdf)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) stop(paste(output, collapse = "\n"))
}

save_ggplot <- function(plot, stem, width, height) {
  final_pdf <- file.path(fig_dir, paste0(stem, ".pdf"))
  raw_pdf <- tempfile(fileext = ".pdf")
  ggsave(raw_pdf, plot, width = width, height = height, device = cm_cairo_pdf)
  convert_pdf_for_latex(raw_pdf, final_pdf)
  unlink(raw_pdf)
}

save_ggplot(p_discrepancy, "fig01_gpd_discrepancy", 8.1, 3.8)
save_ggplot(p_rmse, "fig02_gpd_rmse", 8.1, 4.3)

draw_diagnostics <- function() {
  grid::grid.newpage()
  joint_plot <- ggplotGrob(p_joint)
  location_plot <- ggplotGrob(p_location)
  aligned_widths <- grid::unit.pmax(joint_plot$widths, location_plot$widths)
  aligned_heights <- grid::unit.pmax(joint_plot$heights, location_plot$heights)
  joint_plot$widths <- aligned_widths
  location_plot$widths <- aligned_widths
  joint_plot$heights <- aligned_heights
  location_plot$heights <- aligned_heights
  event_plot <- ggplotGrob(p_events)

  layout <- grid::grid.layout(
    nrow = 2,
    ncol = 2,
    heights = grid::unit(c(0.54, 0.46), "null"),
    widths = grid::unit(c(0.5, 0.5), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(joint_plot)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.draw(location_plot)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  grid::grid.draw(event_plot)
  grid::upViewport(2)
}

diagnostic_raw_pdf <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(
  diagnostic_raw_pdf,
  width = 8.1, height = 7.0,
  family = figure_family, symbolfamily = figure_family
)
draw_diagnostics()
grDevices::dev.off()
convert_pdf_for_latex(
  diagnostic_raw_pdf,
  file.path(fig_dir, "fig03_gpd_diagnostics.pdf")
)
unlink(diagnostic_raw_pdf)

diagnostic_summary <- data.frame(
  n = diagnostic_n,
  growth = diagnostic_growth,
  xi_true = diagnostic$xi_true[1],
  exact_scale_mean = mean(diagnostic$z_sigma_exact),
  exact_shape_mean = mean(diagnostic$z_xi_exact),
  exact_scale_variance = var(diagnostic$z_sigma_exact),
  exact_shape_variance = var(diagnostic$z_xi_exact),
  exact_covariance = cov(diagnostic$z_sigma_exact, diagnostic$z_xi_exact),
  approximate_scale_mean = mean(diagnostic_approximate$z_sigma_approximate),
  approximate_shape_mean = mean(diagnostic_approximate$z_xi_approximate),
  approximate_scale_variance = var(diagnostic_approximate$z_sigma_approximate),
  approximate_shape_variance = var(diagnostic_approximate$z_xi_approximate),
  approximate_covariance = cov(
    diagnostic_approximate$z_sigma_approximate,
    diagnostic_approximate$z_xi_approximate
  ),
  location_mean = mean(diagnostic$location_scaled),
  location_variance = var(diagnostic$location_scaled),
  location_ks_exp1 = suppressWarnings(ks.test(diagnostic$location_scaled, "pexp")$statistic),
  raw_root_interior_rate = mean(diagnostic$approximate_interior == 1),
  pwm_valid_rate = mean(diagnostic$pwm_valid == 1)
)

# Save only summaries that underlie a figure or a numerical statement in the
# manuscript.
write.csv(summary, file.path(out_dir, "gpd_summary.csv"), row.names = FALSE)
write.csv(rmse_summary, file.path(out_dir, "rmse_summary.csv"), row.names = FALSE)
write.csv(event_summary, file.path(out_dir, "event_summary.csv"), row.names = FALSE)
write.csv(diagnostic_summary, file.path(out_dir, "diagnostic_summary.csv"), row.names = FALSE)

config <- c(
  "seed=850000 + 1000*growth_index + n_index",
  paste0("B=", B),
  paste0("chunk_size=", chunk_size),
  paste0("n_grid=", paste(n_grid, collapse = ",")),
  paste0("growth_grid=", paste(growth_grid, collapse = ",")),
  paste0("shape_coefficient=", shape_coefficient),
  paste0("shape_upper_bound=", shape_upper_bound),
  paste0("diagnostic_cell=(n,gamma)=", diagnostic_n, ",", diagnostic_growth),
  "true location=0",
  "true scale=1",
  "benchmark=Hosking-Wallis PWM formulas applied to excesses over the sample minimum",
  "approximate estimator=raw surrogate root when it is in the restricted parameter space; otherwise shifted-exponential boundary estimator (M1,0)",
  "approximate summaries=all replications with a valid exact MLE",
  "PWM RMSE ratios=PWM and exact MLE evaluated on the same PWM-valid replications"
)
writeLines(config, file.path(out_dir, "simulation_config.txt"))
message("Done. Outputs written to ", out_dir, " and ", fig_dir)
