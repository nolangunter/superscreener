#' Simulation design of the SuperScreener chapter
#'
#' Builds a fixed population: covariates are multivariate normal with a
#' latent-factor correlation matrix, and the binary outcome follows
#' \eqn{P(Y=1|W) = \mathrm{expit}\{a + s Z(W)\}}, where \eqn{Z} is the outcome
#' surface \eqn{f} standardized to mean 0 and variance 1. The intercept `a` and
#' slope `s` are calibrated by Monte Carlo so that the prevalence is `prev` and
#' the AUC of the true risk is `auc`.
#'
#' Covariate roles:
#' \describe{
#'   \item{true}{`true1`..`true10` enter \eqn{f}, effects from strong to weak.
#'     Each is its own latent factor; factors have correlation \eqn{0.4^{|k-l|}}.}
#'   \item{proxy}{two per true covariate (`prox_true3_1`, ...), correlated with
#'     it (correlation drawn once from `U(prox_cor)`), outside \eqn{f}.}
#'   \item{noise}{the remaining `p - 30` covariates, in clusters of five with
#'     within-cluster correlation about 0.57.}
#' }
#' A global factor adds 0.15 to every pair. Surfaces: `"nonlinear"` (the base
#' design: linear terms, a steep logistic in `true3`, a step in `true4`, a
#' square in `true5`, and `true1 * true2`), `"linear"`, and
#' `"strong_nonlinear"` (sine, step, products, squares and absolute values).
#'
#' All randomness is drawn from `design_seed`, so a design is the same
#' population on every machine. [chapter_scenarios()] lists the four scenarios
#' of the chapter.
#'
#' @param p Number of covariates (at least 31).
#' @param prox_cor Range of the proxy correlations.
#' @param surface Outcome surface.
#' @param prev Prevalence of the outcome.
#' @param auc AUC of the true risk.
#' @param design_seed Seed for the design.
#' @param mc_n Monte Carlo sample size for the calibration.
#' @param name Optional label.
#' @return A design object (list) for [simulate_data()]. `roles` gives the
#'   role of every covariate and `S0` the true covariates.
#' @examples
#' d <- make_design(p = 50)
#' table(d$roles)
#' @export
make_design <- function(p = 300, prox_cor = c(0.5, 0.85),
                        surface = c("nonlinear", "linear", "strong_nonlinear"),
                        prev = 0.30, auc = 0.80, design_seed = 20261019, mc_n = 2e5, name = "design") {
  surface <- match.arg(surface)
  n_true <- 10L; prox_per_true <- 2L
  stopifnot(p > n_true + prox_per_true * n_true)
  rho_global <- 0.15; rho_between <- 0.4; rho_noise <- 0.5; cluster <- 5
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv))
  set.seed(design_seed)

  tru <- sprintf("true%d", seq_len(n_true))
  n_prox <- prox_per_true * n_true
  n_noise <- p - n_true - n_prox
  anchor <- rep(tru, length.out = n_prox)
  prox <- sprintf("prox_%s_%d", anchor, stats::ave(seq_along(anchor), anchor, FUN = seq_along))
  noise <- sprintf("noise%03d", seq_len(n_noise))
  vars <- c(tru, prox, noise)
  roles <- stats::setNames(c(rep("true", n_true), rep("proxy", n_prox), rep("noise", n_noise)), vars)

  ## latent-factor correlation matrix
  q <- ceiling(n_noise / cluster)
  L <- matrix(0, p, n_true + q, dimnames = list(vars, NULL))
  for (k in seq_len(n_true)) L[tru[k], k] <- 1
  target <- stats::setNames(stats::runif(n_prox, prox_cor[1], prox_cor[2]), prox)
  for (k in seq_len(n_prox))
    L[prox[k], match(anchor[k], tru)] <- min((target[k] - rho_global) / (1 - rho_global), 0.99)
  for (k in seq_len(n_noise)) L[noise[k], n_true + ceiling(k / cluster)] <- sqrt(rho_noise)
  Phi <- diag(n_true + q)
  Phi[seq_len(n_true), seq_len(n_true)] <- rho_between^abs(outer(seq_len(n_true), seq_len(n_true), "-"))
  CU <- L %*% Phi %*% t(L); diag(CU) <- 1
  R <- rho_global + (1 - rho_global) * CU; diag(R) <- 1
  dimnames(R) <- list(vars, vars)

  ## outcome surface
  beta <- stats::setNames(c(1, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1), tru)
  z <- stats::rnorm(1e6)
  h3 <- function(x) stats::plogis(3 * x); h4 <- function(x) as.numeric(x > 0.5)
  c3 <- c(mean(h3(z)), stats::sd(h3(z))); c4 <- c(mean(h4(z)), stats::sd(h4(z)))
  rm(z)
  f <- if (surface == "strong_nonlinear") .strong_nonlinear_f(tru) else function(Ws) {
    Ws <- as.matrix(Ws)
    w <- function(j) Ws[, tru[j]]
    lin <- if (surface == "linear") seq_len(n_true) else setdiff(seq_len(n_true), c(3, 4))
    out <- as.numeric(Ws[, tru[lin], drop = FALSE] %*% beta[lin]) + 0.5 * w(1) * w(2)
    if (surface == "nonlinear")
      out <- out + beta[3] * (h3(w(3)) - c3[1]) / c3[2] + beta[4] * (h4(w(4)) - c4[1]) / c4[2] + 0.4 * (w(5)^2 - 1)
    out
  }

  ## calibrate intercept a and slope s
  set.seed(design_seed + 1)
  Wt <- matrix(stats::rnorm(mc_n * n_true), mc_n) %*% chol(R[tru, tru]); colnames(Wt) <- tru
  fW <- f(Wt); f_mean <- mean(fW); f_sd <- stats::sd(fW)
  fz <- (fW - f_mean) / f_sd
  calib <- function(s) {
    a <- stats::uniroot(function(a) mean(stats::plogis(a + s * fz)) - prev, c(-20, 20))$root
    pr <- stats::plogis(a + s * fz)
    y <- stats::rbinom(mc_n, 1, pr)
    c(a = a, auc = auc(y, pr))
  }
  s <- stats::uniroot(function(s) calib(s)["auc"] - auc, c(0.3, 6), tol = 1e-3)$root
  cal <- calib(s)
  rm(Wt, fW, fz, L, CU, Phi)

  list(name = name, p = p, vars = vars, roles = roles, S0 = tru, anchor = stats::setNames(anchor, prox),
       prox_target_cor = target, R = R, Rc = chol(R), beta = beta, surface = surface, f = f,
       binary = list(a = unname(cal["a"]), s = s, f_mean = f_mean, f_sd = f_sd, prev = prev,
                     auc = unname(cal["auc"])),
       prox_cor = prox_cor, design_seed = design_seed)
}

#' Simulate one dataset from a design
#'
#' @param n Sample size.
#' @param design A design from [make_design()].
#' @param seed Optional seed.
#' @return A list with `X` (data frame of covariates), `Y` (binary outcome)
#'   and `p` (true risk \eqn{P(Y=1|W)} of every row).
#' @export
simulate_data <- function(n, design, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- matrix(stats::rnorm(n * design$p), n) %*% design$Rc
  colnames(W) <- design$vars
  b <- design$binary
  pr <- stats::plogis(b$a + b$s * (design$f(W[, design$S0, drop = FALSE]) - b$f_mean) / b$f_sd)
  list(X = as.data.frame(W), Y = stats::rbinom(n, 1, pr), p = pr)
}

#' The four simulation scenarios of the chapter
#'
#' @return A named list; each element gives `n`, `p`, `prox_cor` and `surface`.
#' \describe{
#'   \item{n500}{moderate sample size: n = 500, p = 300}
#'   \item{prox_high}{high proxy correlation U(0.90, 0.97): n = 1000, p = 300}
#'   \item{small_n}{very small sample: n = 100, p = 2000}
#'   \item{strong_nl}{strongly nonlinear surface: n = 200, p = 500}
#' }
#' @export
chapter_scenarios <- function() {
  list(n500 = list(n = 500, p = 300, prox_cor = c(0.5, 0.85), surface = "nonlinear"),
       prox_high = list(n = 1000, p = 300, prox_cor = c(0.90, 0.97), surface = "nonlinear"),
       small_n = list(n = 100, p = 2000, prox_cor = c(0.5, 0.85), surface = "nonlinear"),
       strong_nl = list(n = 200, p = 500, prox_cor = c(0.5, 0.85), surface = "strong_nonlinear"))
}

#' @noRd
.strong_nonlinear_f <- function(tru) function(Ws) {
  Ws <- as.matrix(Ws); w <- function(j) Ws[, tru[j]]
  1.0 * sin(pi * w(1)) + 0.8 * (w(2) > 0.5) + 0.7 * w(3) * w(4) + 0.6 * (w(5)^2 - 1) +
    0.5 * (abs(w(6)) - sqrt(2 / pi)) + 0.4 * stats::plogis(4 * w(7)) + 0.4 * w(8) * (w(9) > 0) +
    0.3 * cos(2 * w(10))
}
