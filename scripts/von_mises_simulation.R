#!/usr/bin/env Rscript

# Monte Carlo experiments for the von Mises example in the paper.
#
# The script samples from the exact circular model, computes the exact MLE, the
# surrogate-root estimator, and the Abeyasekera--Collett benchmark, and creates
# the three figures used in the manuscript.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

options(scipen = 4)

# Output locations can be overridden without changing the script.
out_dir <- Sys.getenv("VONMISES_SIM_OUT_DIR", "output/summaries/von_mises")
fig_dir <- Sys.getenv("VONMISES_SIM_FIG_DIR", "output/figures/von_mises")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("VONMISES_SIM_B", "2000"))
chunk_size <- as.integer(Sys.getenv("VONMISES_SIM_CHUNK", "100"))
quick <- identical(Sys.getenv("VONMISES_SIM_QUICK", "0"), "1")

if (quick) {
  n_grid <- c(250L, 1000L)
} else {
  n_grid <- c(250L, 500L, 1000L, 2500L, 5000L, 10000L, 25000L)
}
growth_grid <- c(0.15, 0.25, 0.40)
diagnostic_n <- max(n_grid)
diagnostic_growth <- 0.40
figure_family <- Sys.getenv("FIGURE_FONT", "CMU Serif")

blue <- "#0072B2"
orange <- "#E69F00"
purple <- "#CC79A7"
green <- "#009E73"
grey <- "#4D4D4D"
growth_colours <- c(
  "0.15" = blue,
  "0.25" = orange,
  "0.4" = purple
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

wrap_angle <- function(x) (x + pi) %% (2 * pi) - pi

# Best--Fisher rejection sampler, vectorized in batches. Angles are returned on
# [0, 2*pi), matching the circular sufficient statistics used below.
r_vonmises <- function(n, kappa, mu = 0) {
  if (!is.finite(kappa) || kappa <= 0) stop("kappa must be positive")
  r <- 1 + sqrt(1 + 4 * kappa^2)
  rho <- (r - sqrt(2 * r)) / (2 * kappa)
  s <- (1 + rho^2) / (2 * rho)
  answer <- numeric(n)
  filled <- 0L
  while (filled < n) {
    m <- max(32L, ceiling(1.25 * (n - filled)))
    u1 <- runif(m)
    u2 <- runif(m)
    u3 <- runif(m)
    z <- cos(pi * u1)
    f <- (1 + s * z) / (s + z)
    c_value <- kappa * (s - f)
    accept <- (u2 < c_value * (2 - c_value)) |
      (log(c_value / u2) + 1 - c_value >= 0)
    values <- ifelse(
      u3[accept] > 0.5,
      acos(f[accept]),
      -acos(f[accept])
    ) + mu
    take <- min(length(values), n - filled)
    if (take > 0L) {
      answer[(filled + 1L):(filled + take)] <- values[seq_len(take)]
    }
    filled <- filled + take
  }
  answer %% (2 * pi)
}

bessel_ratio <- function(kappa) {
  besselI(kappa, nu = 1, expon.scaled = TRUE) /
    besselI(kappa, nu = 0, expon.scaled = TRUE)
}

# Vectorized Newton inversion of A(kappa)=I_1(kappa)/I_0(kappa). The update is
# safeguarded to keep every concentration parameter positive.
inverse_bessel_ratio <- function(R) {
  kappa <- (1 + sqrt(pmax(0, 3 - 2 * R))) / (4 * (1 - R))
  for (iteration in seq_len(12L)) {
    a <- bessel_ratio(kappa)
    derivative <- 1 - a / kappa - a^2
    candidate <- kappa - (a - R) / derivative
    bad <- !is.finite(candidate) | candidate <= 0
    candidate[bad] <- kappa[bad] / 2
    kappa <- candidate
  }
  kappa
}

rmse <- function(x, truth) sqrt(mean((x - truth)^2, na.rm = TRUE))

# Simulate and estimate one (n, gamma) cell in memory-bounded chunks.
simulate_cell <- function(n, growth, B, chunk_size, seed) {
  set.seed(seed)
  kappa <- n^growth
  pieces <- vector("list", ceiling(B / chunk_size))
  piece_index <- 1L

  for (first in seq.int(1L, B, by = chunk_size)) {
    b <- min(chunk_size, B - first + 1L)
    x <- matrix(
      r_vonmises(b * n, kappa = kappa, mu = 0),
      nrow = b,
      ncol = n
    )
    cbar <- rowMeans(cos(x))
    sbar <- rowMeans(sin(x))
    resultant <- sqrt(cbar^2 + sbar^2)
    mu_hat <- atan2(sbar, cbar)

    kappa_hat <- inverse_bessel_ratio(resultant)
    kappa_tilde <-
      (1 + sqrt(pmax(0, 3 - 2 * resultant))) /
      (4 * (1 - resultant))

    centered <- matrix(
      wrap_angle(x - rep(mu_hat, n)),
      nrow = b,
      ncol = n
    )
    centered_mean <- rowMeans(centered)
    centered <- centered - centered_mean
    sum_squares <- rowSums(centered^2)
    kappa_linear <- (n - 3) / sum_squares

    asymptotic_se <- sqrt(2) * kappa / sqrt(n)
    z_mu <- sqrt(n * kappa) * wrap_angle(mu_hat)
    z_kappa_hat <- sqrt(n / 2) * (kappa_hat - kappa) / kappa
    z_kappa_tilde <- sqrt(n / 2) * (kappa_tilde - kappa) / kappa
    z_kappa_linear <- sqrt(n / 2) * (kappa_linear - kappa) / kappa

    pieces[[piece_index]] <- data.frame(
      growth = growth,
      n = n,
      kappa_true = kappa,
      resultant = resultant,
      mu_hat = mu_hat,
      kappa_hat = kappa_hat,
      kappa_tilde = kappa_tilde,
      kappa_linear = kappa_linear,
      score_residual = bessel_ratio(kappa_hat) - resultant,
      discrepancy_se = abs(kappa_tilde - kappa_hat) / asymptotic_se,
      z_mu = z_mu,
      z_kappa_hat = z_kappa_hat,
      z_kappa_tilde = z_kappa_tilde,
      z_kappa_linear = z_kappa_linear
    )
    piece_index <- piece_index + 1L
  }
  bind_rows(pieces)
}

summarize_cell <- function(dat) {
  truth <- dat$kappa_true[1]
  exact_rmse <- rmse(dat$kappa_hat, truth)
  discrepancy_quantiles <- quantile(
    dat$discrepancy_se,
    probs = c(0.05, 0.50, 0.95),
    na.rm = TRUE
  )
  data.frame(
    growth = dat$growth[1],
    n = dat$n[1],
    kappa_true = truth,
    discrepancy_q05 = unname(discrepancy_quantiles[1]),
    discrepancy_median = unname(discrepancy_quantiles[2]),
    discrepancy_q95 = unname(discrepancy_quantiles[3]),
    exact_bias = mean(dat$kappa_hat - truth, na.rm = TRUE),
    exact_rmse = exact_rmse,
    approximate_bias = mean(dat$kappa_tilde - truth, na.rm = TRUE),
    approximate_rmse = rmse(dat$kappa_tilde, truth),
    approximate_rmse_ratio = rmse(dat$kappa_tilde, truth) / exact_rmse,
    linear_bias = mean(dat$kappa_linear - truth, na.rm = TRUE),
    linear_rmse = rmse(dat$kappa_linear, truth),
    linear_rmse_ratio = rmse(dat$kappa_linear, truth) / exact_rmse,
    max_abs_score_residual = max(abs(dat$score_residual), na.rm = TRUE)
  )
}

# Deterministic cell-specific seeds make the results reproducible regardless of
# the order in which the grid is evaluated.
message(
  "Von Mises study: ", length(n_grid) * length(growth_grid),
  " cells; B = ", B
)

cells <- vector("list", length(n_grid) * length(growth_grid))
cell <- 1L
for (growth_index in seq_along(growth_grid)) {
  growth <- growth_grid[growth_index]
  for (n_index in seq_along(n_grid)) {
    n <- n_grid[n_index]
    message("  gamma = ", growth, ", n = ", n)
    cells[[cell]] <- simulate_cell(
      n, growth, B, chunk_size,
      seed = 830000L + 1000L * growth_index + n_index
    )
    cell <- cell + 1L
  }
}

replications <- bind_rows(cells)
summary <- bind_rows(lapply(cells, summarize_cell))
summary <- summary[order(summary$growth, summary$n), ]
summary$growth_factor <- factor(summary$growth)

rmse_summary <- bind_rows(
  summary %>% transmute(
    growth, growth_factor, n,
    method = "Approximate MLE",
    ratio = approximate_rmse_ratio
  ),
  summary %>% transmute(
    growth, growth_factor, n,
    method = "Abeyasekera-Collett",
    ratio = linear_rmse_ratio
  )
)
rmse_summary$method <- factor(
  rmse_summary$method,
  levels = c("Approximate MLE", "Abeyasekera-Collett")
)

p_discrepancy <- ggplot(
  summary,
  aes(x = n, y = discrepancy_median, colour = growth_factor, group = growth_factor)
) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.35, linetype = "longdash") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.25) +
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(
    breaks = 10^seq(-4, 2, by = 1),
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
  "Abeyasekera-Collett" = "longdash"
)
method_shapes <- c(
  "Approximate MLE" = 16,
  "Abeyasekera-Collett" = 17
)

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
  scale_x_log10(
    breaks = c(250, 1000, 5000, 25000),
    labels = c("250", "1,000", "5,000", "25,000")
  ) +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20),
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

density_data <- bind_rows(
  data.frame(estimator = "Exact MLE", z = diagnostic$z_kappa_hat),
  data.frame(estimator = "Approximate MLE", z = diagnostic$z_kappa_tilde)
)

p_density <- ggplot(density_data, aes(x = z, colour = estimator, linetype = estimator)) +
  geom_density(linewidth = 0.75, adjust = 1.1, key_glyph = "path") +
  stat_function(
    fun = dnorm,
    aes(colour = "Standard normal", linetype = "Standard normal"),
    linewidth = 0.65,
    key_glyph = "path"
  ) +
  scale_colour_manual(
    values = c(
      "Exact MLE" = blue,
      "Approximate MLE" = orange,
      "Standard normal" = grey
    ),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "Exact MLE" = "solid",
      "Approximate MLE" = "dashed",
      "Standard normal" = "dotted"
    ),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-4, 4)) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x = "Standardized concentration error",
    y = "Density",
    title = "(a) Von Mises: concentration limit"
  ) +
  paper_theme +
  theme(plot.title = element_text(size = 11, face = "plain"))

circle_probability <- 0.90
circle_radius <- sqrt(qchisq(circle_probability, df = 2))
angles <- seq(0, 2 * pi, length.out = 300)
circle <- data.frame(
  x = circle_radius * cos(angles),
  y = circle_radius * sin(angles)
)

p_joint <- ggplot(diagnostic, aes(x = z_mu, y = z_kappa_tilde)) +
  geom_point(colour = blue, alpha = 0.20, size = 0.65) +
  geom_path(data = circle, aes(x = x, y = y), inherit.aes = FALSE,
            colour = grey, linetype = "dashed", linewidth = 0.6) +
  coord_equal(xlim = c(-4, 4), ylim = c(-4, 4)) +
  labs(
    x = expression(sqrt(italic(n) * kappa[n]) * (hat(mu)[n] - mu)),
    y = expression(sqrt(italic(n) / 2) * (tilde(kappa)[n] - kappa[n]) / kappa[n]),
    title = "(b) Von Mises: joint limit"
  ) +
  paper_theme +
  theme(legend.position = "none", plot.title = element_text(size = 11, face = "plain"))

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

save_ggplot(p_discrepancy, "fig01_vonmises_discrepancy", 8.1, 3.8)
save_ggplot(p_rmse, "fig02_vonmises_rmse", 8.1, 4.3)

draw_diagnostics <- function() {
  grid::grid.newpage()
  density_table <- ggplotGrob(p_density)
  legend_index <- which(vapply(
    density_table$grobs,
    function(grob) identical(grob$name, "guide-box"),
    logical(1)
  ))
  shared_legend <- density_table$grobs[[legend_index]]

  density_plot <- ggplotGrob(p_density + theme(legend.position = "none"))
  joint_plot <- ggplotGrob(p_joint)
  aligned_widths <- grid::unit.pmax(density_plot$widths, joint_plot$widths)
  aligned_heights <- grid::unit.pmax(density_plot$heights, joint_plot$heights)
  density_plot$widths <- aligned_widths
  joint_plot$widths <- aligned_widths
  density_plot$heights <- aligned_heights
  joint_plot$heights <- aligned_heights

  layout <- grid::grid.layout(
    nrow = 2,
    ncol = 2,
    heights = grid::unit(c(0.84, 0.16), "null"),
    widths = grid::unit(c(0.5, 0.5), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(density_plot)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.draw(joint_plot)
  grid::upViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  grid::grid.draw(shared_legend)
  grid::upViewport(2)
}

diagnostic_raw_pdf <- tempfile(fileext = ".pdf")
grDevices::cairo_pdf(
  diagnostic_raw_pdf,
  width = 8.1, height = 4.0,
  family = figure_family, symbolfamily = figure_family
)
draw_diagnostics()
grDevices::dev.off()
convert_pdf_for_latex(
  diagnostic_raw_pdf,
  file.path(fig_dir, "fig03_vonmises_diagnostics.pdf")
)
unlink(diagnostic_raw_pdf)

diagnostic_summary <- data.frame(
  n = diagnostic_n,
  growth = diagnostic_growth,
  kappa = diagnostic$kappa_true[1],
  exact_concentration_mean = mean(diagnostic$z_kappa_hat),
  exact_concentration_variance = var(diagnostic$z_kappa_hat),
  approximate_concentration_mean = mean(diagnostic$z_kappa_tilde),
  approximate_concentration_variance = var(diagnostic$z_kappa_tilde),
  location_mean = mean(diagnostic$z_mu),
  location_variance = var(diagnostic$z_mu),
  approximate_joint_correlation = cor(diagnostic$z_mu, diagnostic$z_kappa_tilde)
)

# Save only summaries that underlie a figure or a numerical statement in the
# manuscript.
write.csv(summary, file.path(out_dir, "vonmises_summary.csv"), row.names = FALSE)
write.csv(rmse_summary, file.path(out_dir, "rmse_summary.csv"), row.names = FALSE)
write.csv(diagnostic_summary, file.path(out_dir, "diagnostic_summary.csv"), row.names = FALSE)

config <- c(
  "seed=830000 + 1000*growth_index + n_index",
  paste0("B=", B),
  paste0("chunk_size=", chunk_size),
  paste0("n_grid=", paste(n_grid, collapse = ",")),
  paste0("growth_grid=", paste(growth_grid, collapse = ",")),
  paste0("diagnostic_cell=(n,gamma)=", diagnostic_n, ",", diagnostic_growth),
  "true location=0",
  "exact concentration MLE=vectorized Newton inversion of I1/I0",
  "benchmark=Abeyasekera-Collett linear estimator after circular unwrapping"
)
writeLines(config, file.path(out_dir, "simulation_config.txt"))
message("Done. Outputs written to ", out_dir, " and ", fig_dir)
