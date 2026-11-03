#!/usr/bin/env Rscript

# Black cherry tree application for the gamma-regression example.
#
# The data are distributed with R. The script fits the conventional gamma GLM,
# computes the raw and bias-corrected explicit estimators, and prints the
# estimates, Wald intervals, fitted-mean comparisons, and small-dispersion
# diagnostics reported in the paper.

# Fit the mean model with a gamma response and logarithmic link.
data(trees, package = "datasets")

fit <- glm(
  Volume ~ Girth + Height,
  data = trees,
  family = Gamma(link = "log"),
  control = glm.control(epsilon = 1e-12, maxit = 100)
)

y <- model.response(model.frame(fit))
X <- model.matrix(fit)
n <- nrow(X)
k <- ncol(X)

# The explicit coefficient estimator is least squares on the log responses.
# The dispersion correction and intercept correction use the residual degrees
# of freedom and the second-order small-dispersion expansion.
beta_raw <- drop(qr.solve(X, log(y)))
log_residual <- drop(log(y) - X %*% beta_raw)
phi_raw <- mean(log_residual^2)
v_check <- phi_raw / (1 - k / n)
phi_bc <- v_check - 0.5 * v_check^2
beta_bc <- beta_raw
beta_bc[1] <- beta_bc[1] + 0.5 * v_check - v_check^2 / 6

# Solve the full-likelihood dispersion equation at the GLM coefficient root.
z_hat <- y / fitted(fit)
target <- -mean(log(z_hat))
g <- function(phi) -log(phi) - digamma(1 / phi) - target
phi_mle <- uniroot(g, interval = c(1e-10, 1e3), tol = 1e-14)$root

phi_pearson <- summary(fit)$dispersion
phi_deviance <- deviance(fit) / df.residual(fit)

# Assemble the coefficient and dispersion comparisons used in the table.
comparison <- data.frame(
  parameter = c(names(coef(fit)), "Dispersion"),
  glm = c(unname(coef(fit)), phi_pearson),
  raw_explicit = c(unname(beta_raw), phi_raw),
  bias_corrected = c(unname(beta_bc), phi_bc)
)
comparison$raw_difference <- comparison$raw_explicit - comparison$glm
comparison$corrected_difference <- comparison$bias_corrected - comparison$glm
comparison$raw_difference_over_glm_se <- c(
  comparison$raw_difference[seq_len(k)] / coef(summary(fit))[, "Std. Error"],
  NA_real_
)
comparison$corrected_difference_over_glm_se <- c(
  comparison$corrected_difference[seq_len(k)] / coef(summary(fit))[, "Std. Error"],
  NA_real_
)

mu_glm <- fitted(fit)
mu_raw <- drop(exp(X %*% beta_raw))
mu_bc <- drop(exp(X %*% beta_bc))

# First-order covariance matrices and 95% Wald intervals.
vcov_raw <- phi_raw * chol2inv(chol(crossprod(X)))
vcov_bc <- phi_bc * chol2inv(chol(crossprod(X)))
se_raw <- sqrt(diag(vcov_raw))
se_bc <- sqrt(diag(vcov_bc))
se_glm <- coef(summary(fit))[, "Std. Error"]
z_975 <- qnorm(0.975)
ci_comparison <- data.frame(
  parameter = names(coef(fit)),
  glm_estimate = unname(coef(fit)),
  glm_se = unname(se_glm),
  glm_lower = unname(coef(fit) - z_975 * se_glm),
  glm_upper = unname(coef(fit) + z_975 * se_glm),
  raw_estimate = unname(beta_raw),
  raw_se = unname(se_raw),
  raw_lower = unname(beta_raw - z_975 * se_raw),
  raw_upper = unname(beta_raw + z_975 * se_raw),
  bc_estimate = unname(beta_bc),
  bc_se = unname(se_bc),
  bc_lower = unname(beta_bc - z_975 * se_bc),
  bc_upper = unname(beta_bc + z_975 * se_bc)
)

se_phi_glm <- sqrt(2 / n) * phi_pearson
se_phi_raw <- sqrt(2 / n) * phi_raw
se_phi_bc <- sqrt(2 / n) * phi_bc
se_phi_mle <- phi_mle^2 /
  sqrt(n * (trigamma(1 / phi_mle) - phi_mle))

# Compare the exact gamma deviance contribution with its quadratic expansion.
h <- function(z) z - 1 - log(z)
exact_deviance_quantity <- mean(2 * h(exp(log_residual)))
relative_deviance_difference <-
  exact_deviance_quantity / phi_raw - 1

# Print every numerical quantity cited in the application.
cat("n =", n, "; k =", k, "\n")
cat("GLM converged:", fit$converged, "; iterations:", fit$iter, "\n\n")
print(comparison, row.names = FALSE, digits = 10)
cat("\nDispersion estimates\n")
print(
  data.frame(
    method = c(
      "Full likelihood",
      "Raw explicit",
      "Bias-corrected explicit",
      "Pearson / residual df",
      "Deviance / residual df"
    ),
    estimate = c(phi_mle, phi_raw, phi_bc, phi_pearson, phi_deviance)
  ),
  row.names = FALSE,
  digits = 10
)
cat("\nFitted-mean comparison\n")
print(
  data.frame(
    method = c("Raw explicit", "Bias-corrected explicit"),
    mean_absolute_relative_difference = c(
      mean(abs(mu_raw / mu_glm - 1)),
      mean(abs(mu_bc / mu_glm - 1))
    ),
    maximum_absolute_relative_difference = c(
      max(abs(mu_raw / mu_glm - 1)),
      max(abs(mu_bc / mu_glm - 1))
    ),
    correlation = c(
      cor(mu_raw, mu_glm),
      cor(mu_bc, mu_glm)
    )
  ),
  row.names = FALSE,
  digits = 10
)
cat("\nCoefficient standard errors and 95% Wald intervals\n")
print(ci_comparison, row.names = FALSE, digits = 10)
cat("\nDispersion standard errors\n")
print(
  data.frame(
    method = c(
      "Full likelihood Fisher information",
      "GLM Pearson first-order",
      "Raw explicit first-order",
      "Bias-corrected first-order"
    ),
    estimate = c(phi_mle, phi_pearson, phi_raw, phi_bc),
    standard_error = c(se_phi_mle, se_phi_glm, se_phi_raw, se_phi_bc),
    lower = pmax(
      0,
      c(phi_mle, phi_pearson, phi_raw, phi_bc) -
        z_975 * c(se_phi_mle, se_phi_glm, se_phi_raw, se_phi_bc)
    ),
    upper = c(phi_mle, phi_pearson, phi_raw, phi_bc) +
      z_975 * c(se_phi_mle, se_phi_glm, se_phi_raw, se_phi_bc)
  ),
  row.names = FALSE,
  digits = 10
)
cat("\nSmall-dispersion diagnostics from the raw explicit fit\n")
print(
  c(
    raw_dispersion = phi_raw,
    raw_dispersion_sqrt = sqrt(phi_raw),
    maximum_absolute_log_residual = max(abs(log_residual)),
    exact_deviance_quantity = exact_deviance_quantity,
    relative_deviance_difference = relative_deviance_difference
  ),
  digits = 10
)
