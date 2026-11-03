#!/usr/bin/env Rscript

# Monte Carlo experiments for the Zeta-distribution example in the paper.
#
# The script simulates the exact discrete model, computes the exact and
# surrogate-root estimators, and produces the three Zeta figures together with
# the compact summaries used in the manuscript.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

options(scipen = 4)

# Output locations can be overridden without changing the script.
out_dir <- Sys.getenv("ZETA_SIM_OUT_DIR", "output/summaries/zeta")
fig_dir <- Sys.getenv("ZETA_SIM_FIG_DIR", "output/figures/zeta")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Full runs use the manuscript grid and 2,000 replications per cell. Quick mode
# retains all three regimes on a two-point sample-size grid.
B <- as.integer(Sys.getenv("ZETA_SIM_B", "2000"))
chunk_size <- as.integer(Sys.getenv("ZETA_SIM_CHUNK", "100"))
quick <- identical(Sys.getenv("ZETA_SIM_QUICK", "0"), "1")
n_grid <- if (quick) {
  c(250L, 1000L)
} else {
  c(250L, 500L, 1000L, 2500L, 5000L, 10000L, 25000L)
}
growth_grid <- c(0.40, 0.50, 0.70)
diagnostic_n <- max(n_grid)
diagnostic_growth <- 0.70
log2_const <- log(2)
figure_family <- Sys.getenv("FIGURE_FONT", "CMU Serif")

blue <- "#0072B2"
orange <- "#D55E00"
green <- "#009E73"
purple <- "#CC79A7"
grey <- "#4D4D4D"
growth_colours <- c(
  "0.4" = blue,
  "0.5" = "#E69F00",
  "0.7" = purple
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

# Euler--Maclaurin evaluation of zeta(s) and -zeta'(s), including terms
# through B_8. With 1,000 direct terms, the numerical error is negligible
# relative to the Monte Carlo error over the simulation grid.
zeta_moments <- function(s, n_terms = 1000L) {
  stopifnot(length(s) == 1L, is.finite(s), s > 1)
  k <- seq_len(n_terms - 1L)
  log_k <- log(k)
  weights <- exp(-s * log_k)
  z <- sum(weights)
  minus_z_prime <- sum(log_k * weights)
  N <- as.numeric(n_terms)
  log_N <- log(N)

  integral <- N^(1 - s) / (s - 1)
  z <- z + integral
  minus_z_prime <- minus_z_prime + integral * (log_N + 1 / (s - 1))

  endpoint <- 0.5 * N^(-s)
  z <- z + endpoint
  minus_z_prime <- minus_z_prime + endpoint * log_N

  add_bernoulli_term <- function(coefficient, order, power) {
    factors <- s + seq.int(0, order - 1L)
    term <- coefficient * prod(factors) * N^(-s - power)
    h <- sum(1 / factors)
    c(zeta = term, minus_derivative = term * (log_N - h))
  }
  for (term in list(
      add_bernoulli_term(1 / 12, 1L, 1L),
      add_bernoulli_term(-1 / 720, 3L, 3L),
      add_bernoulli_term(1 / 30240, 5L, 5L),
      add_bernoulli_term(-1 / 1209600, 7L, 7L))) {
    z <- z + term[["zeta"]]
    minus_z_prime <- minus_z_prime + term[["minus_derivative"]]
  }
  c(zeta = z, mean_log = minus_z_prime / z)
}

# The exact MLE solves E_s(log X) = log(G_n). A dense monotone lookup avoids
# evaluating the Dirichlet series separately for every Monte Carlo sample.
build_score_inverse <- function() {
  s_grid <- seq(1.0005, 50, length.out = 30000L)
  mean_log <- vapply(s_grid, function(s) zeta_moments(s)[["mean_log"]], 0.0)
  inverse <- splinefun(rev(log(mean_log)), rev(s_grid), method = "monoH.FC")
  list(
    evaluate = function(y) {
      if (any(!is.finite(y) | y <= 0)) stop("score inversion requires y > 0")
      log_y <- log(y)
      if (any(log_y < min(log(mean_log)) | log_y > max(log(mean_log)))) {
        stop("score value outside the precomputed inversion range")
      }
      inverse(log_y)
    },
    s_grid = s_grid,
    mean_log = mean_log
  )
}

# Rejection sampler for the Zeta(s) distribution conditional on X >= 2.
# Sampling the number of non-ones first preserves the exact iid law without
# allocating a B by n matrix dominated by ones.
rzeta_conditional_ge2 <- function(count, s) {
  if (count == 0L) return(numeric())
  answer <- numeric(count)
  pending <- seq_len(count)
  q_two <- 2^(-s) / (2^(1 - s) - 3^(1 - s))
  while (length(pending) > 0L) {
    y <- 2 * runif(length(pending))^(-1 / (s - 1))
    k <- pmax(2, floor(y))
    q_k <- k^(-s) / (k^(1 - s) - (k + 1)^(1 - s))
    accepted <- runif(length(pending)) <= q_k / q_two
    answer[pending[accepted]] <- k[accepted]
    pending <- pending[!accepted]
  }
  answer
}

simulate_cell <- function(n, growth, B, inverse_score, chunk_size, seed) {
  set.seed(seed)
  u_true <- n^(-growth)
  s_true <- -log(u_true) / log2_const
  zeta_value <- zeta_moments(s_true)[["zeta"]]
  p_nonone <- 1 - 1 / zeta_value
  log_g <- numeric(B)

  starts <- seq.int(1L, B, by = chunk_size)
  for (start in starts) {
    stop_at <- min(B, start + chunk_size - 1L)
    indices <- start:stop_at
    nonone_counts <- rbinom(length(indices), n, p_nonone)
    total_nonones <- sum(nonone_counts)
    if (total_nonones > 0L) {
      values <- rzeta_conditional_ge2(total_nonones, s_true)
      groups <- rep.int(seq_along(indices), nonone_counts)
      grouped_sums <- rowsum(log(values), group = groups, reorder = FALSE)
      log_sums <- numeric(length(indices))
      log_sums[as.integer(rownames(grouped_sums))] <- grouped_sums[, 1]
      log_g[indices] <- log_sums / n
    }
  }

  finite <- log_g > 0
  u_tilde <- s_tilde <- u_hat <- s_hat <- rep(NA_real_, B)
  u_tilde[finite] <- log_g[finite] / log2_const
  s_tilde[finite] <- -log(u_tilde[finite]) / log2_const
  s_hat[finite] <- inverse_score(log_g[finite])
  u_hat[finite] <- 2^(-s_hat[finite])

  score_residual <- rep(NA_real_, B)
  score_residual[finite] <- vapply(s_hat[finite], function(s) {
    zeta_moments(s)[["mean_log"]]
  }, 0.0) - log_g[finite]

  data.frame(
    growth = growth,
    n = n,
    u_true = u_true,
    s_true = s_true,
    log_g = log_g,
    finite_mle = finite,
    approximate_in_parameter_space = finite & s_tilde > 1,
    u_tilde = u_tilde,
    u_hat = u_hat,
    s_tilde = s_tilde,
    s_hat = s_hat,
    score_residual = score_residual,
    stringsAsFactors = FALSE
  )
}

rmse <- function(x, truth) sqrt(mean((x - truth)^2))

summarize_cell <- function(dat) {
  d <- dat[dat$finite_mle, ]
  se_s <- 1 / (log2_const * sqrt(d$n[1] * d$u_true[1]))
  discrepancy_s <- abs(d$s_tilde - d$s_hat) / se_s
  exact_rmse <- rmse(d$s_hat, d$s_true[1])

  data.frame(
    growth = dat$growth[1],
    n = dat$n[1],
    u_true = dat$u_true[1],
    s_true = dat$s_true[1],
    finite_mle_rate = mean(dat$finite_mle),
    approximate_in_parameter_space_rate =
      mean(dat$approximate_in_parameter_space[dat$finite_mle]),
    discrepancy_s_median = median(discrepancy_s),
    discrepancy_s_q05 = unname(quantile(discrepancy_s, 0.05)),
    discrepancy_s_q95 = unname(quantile(discrepancy_s, 0.95)),
    rmse_s_tilde = rmse(d$s_tilde, d$s_true[1]),
    rmse_s_hat = exact_rmse,
    rmse_ratio_tilde = rmse(d$s_tilde, d$s_true[1]) / exact_rmse,
    boundary_rate = mean(!dat$finite_mle),
    boundary_theory = zeta_moments(dat$s_true[1])[["zeta"]]^(-dat$n[1]),
    max_abs_score_residual = max(abs(d$score_residual))
  )
}

lookup <- build_score_inverse()
# Deterministic cell-specific seeds make the simulation independent of the
# order in which the grid points are evaluated.
all_cells <- vector("list", length(growth_grid) * length(n_grid))
cell <- 1L
for (g_index in seq_along(growth_grid)) {
  growth <- growth_grid[g_index]
  for (n_index in seq_along(n_grid)) {
    n <- n_grid[n_index]
    message("Zeta simulation: gamma = ", growth, ", n = ", n, ", B = ", B)
    all_cells[[cell]] <- simulate_cell(
      n, growth, B, lookup$evaluate, chunk_size,
      seed = 910000L + 1000L * g_index + n_index
    )
    cell <- cell + 1L
  }
}

replications <- bind_rows(all_cells)
summary <- bind_rows(lapply(all_cells, summarize_cell))
summary <- summary[order(summary$growth, summary$n), ]

discrepancy_summary <- summary %>% transmute(
  growth, n,
  median = discrepancy_s_median, q05 = discrepancy_s_q05, q95 = discrepancy_s_q95
)
discrepancy_summary$growth_factor <- factor(discrepancy_summary$growth)

rmse_summary <- summary %>% transmute(
  growth, n, ratio = rmse_ratio_tilde
)
rmse_summary$growth_factor <- factor(rmse_summary$growth)

p_discrepancy <- ggplot(
  discrepancy_summary,
  aes(x = n, y = median, colour = growth_factor, group = growth_factor)
) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.35, linetype = "longdash") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.25) +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(
    breaks = 10^seq(-4, 3, by = 1),
    labels = power_of_ten_labels
  ) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  guides(colour = guide_legend(nrow = 1)) +
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
    group = growth_factor
  )
) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.4, linetype = "longdash") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.25, stroke = 0.45) +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(breaks = c(1, 2, 3, 4, 5), labels = scales::label_number(accuracy = 1)) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(
    x = expression(italic(n)),
    y = "RMSE relative to the exact MLE"
  ) +
  paper_theme

diagnostic <- replications[
  replications$n == diagnostic_n &
    abs(replications$growth - diagnostic_growth) < 1e-12 &
    replications$finite_mle,
]

make_qq_data <- function(exact, approximate, scale) {
  probs <- ppoints(length(exact))
  bind_rows(
    data.frame(
      scale = scale,
      estimator = "Exact MLE",
      theoretical = qnorm(probs),
      empirical = sort(exact)
    ),
    data.frame(
      scale = scale,
      estimator = "Approximate MLE",
      theoretical = qnorm(probs),
      empirical = sort(approximate)
    )
  )
}

z_s_hat <- with(diagnostic, log2_const * sqrt(n * u_true) * (s_hat - s_true))
z_s_tilde <- with(diagnostic, log2_const * sqrt(n * u_true) * (s_tilde - s_true))

qq_original <- make_qq_data(z_s_hat, z_s_tilde, "Original parameter scale")

qq_colours <- c("Exact MLE" = blue, "Approximate MLE" = orange)
qq_linetypes <- c("Exact MLE" = "solid", "Approximate MLE" = "dashed")

make_qq_plot <- function(dat, title) {
  ggplot(dat, aes(theoretical, empirical, colour = estimator, linetype = estimator)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(values = qq_colours, name = NULL) +
    scale_linetype_manual(values = qq_linetypes, name = NULL) +
    coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5)) +
    labs(
      x = "Standard normal quantile",
      y = "Empirical quantile",
      title = title
    ) +
    paper_theme +
    theme(plot.title = element_text(size = 11, face = "plain"))
}

p_original_qq <- make_qq_plot(
  qq_original,
  "(a) Original parameter scale: standardized limit"
)

boundary_summary <- summary
boundary_summary$growth_factor <- factor(boundary_summary$growth)
boundary_long <- bind_rows(
  boundary_summary %>% transmute(
    growth_factor, n, source = "Empirical", rate = boundary_rate
  ),
  boundary_summary %>% transmute(
    growth_factor, n, source = "Exact probability", rate = boundary_theory
  )
)
boundary_long$source <- factor(
  boundary_long$source,
  levels = c("Empirical", "Exact probability")
)

p_boundary <- ggplot(
  boundary_long,
  aes(
    x = n,
    y = rate,
    colour = growth_factor,
    linetype = source,
    group = interaction(growth_factor, source)
  )
) +
  geom_line(linewidth = 0.65) +
  geom_point(
    data = boundary_long[boundary_long$source == "Empirical", ],
    aes(x = n, y = rate, colour = growth_factor),
    inherit.aes = FALSE,
    size = 1.25
  ) +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_sqrt(labels = scales::label_percent(accuracy = 0.01)) +
  scale_colour_manual(values = growth_colours, name = expression(gamma)) +
  scale_linetype_manual(
    values = c("Empirical" = "solid", "Exact probability" = "longdash"),
    name = NULL
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 1)
  ) +
  labs(
    x = expression(italic(n)),
    y = "Probability of no finite MLE",
    title = "(b) Boundary rate"
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
  gs <- Sys.which("gs")
  if (!nzchar(gs)) stop("Ghostscript is required to create the final figure PDFs.")
  output <- system2(
    gs,
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

save_ggplot(p_discrepancy, "fig01_zeta_discrepancy", 6.2, 4.0)
save_ggplot(p_rmse, "fig02_zeta_rmse", 6.2, 4.0)

draw_diagnostics <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(0.5, 0.5), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(p_original_qq, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_boundary, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::upViewport()
}

diagnostic_raw_pdf <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(
  diagnostic_raw_pdf,
  width = 7.2, height = 7.2,
  family = figure_family, symbolfamily = figure_family
)
draw_diagnostics()
grDevices::dev.off()
convert_pdf_for_latex(
  diagnostic_raw_pdf,
  file.path(fig_dir, "fig03_zeta_diagnostics.pdf")
)
unlink(diagnostic_raw_pdf)

# Save the cell summaries that underlie the figures and reported diagnostics.
write.csv(summary, file.path(out_dir, "zeta_summary.csv"), row.names = FALSE)
write.csv(discrepancy_summary, file.path(out_dir, "discrepancy_summary.csv"), row.names = FALSE)
write.csv(rmse_summary, file.path(out_dir, "rmse_summary.csv"), row.names = FALSE)

config <- c(
  "seed=910000 + 1000*growth_index + n_index",
  paste0("B=", B),
  paste0("chunk_size=", chunk_size),
  paste0("n_grid=", paste(n_grid, collapse = ",")),
  paste0("growth_grid=", paste(growth_grid, collapse = ",")),
  paste0("diagnostic_cell=(n,gamma)=", diagnostic_n, ",", diagnostic_growth),
  "exact root=monotone interpolation of an Euler--Maclaurin Dirichlet-series score",
  "sampler=exact conditional Pareto-envelope rejection sampler on {2,3,...}"
)
writeLines(config, file.path(out_dir, "simulation_config.txt"))
message("Done. Outputs written to ", out_dir, " and ", fig_dir)
