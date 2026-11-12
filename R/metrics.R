#' Prediction metrics for a binary outcome
#'
#' `brier()` is the mean squared error of the predicted probabilities (the
#' selection loss of SuperScreener). `neg_loglik()` is the mean negative
#' log-likelihood with predictions truncated to `[1e-6, 1 - 1e-6]`. `auc()` is
#' the area under the ROC curve.
#' @param y Binary outcome (0/1).
#' @param p Predicted probabilities.
#' @name metrics
NULL

#' @rdname metrics
#' @export
brier <- function(y, p) mean((y - p)^2)

#' @rdname metrics
#' @export
neg_loglik <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

#' @rdname metrics
#' @export
auc <- function(y, p) {
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- length(y) - n1
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Cross-validated AUC with an influence-curve confidence interval
#'
#' The cross-validated AUC is the mean of the AUCs of the validation folds. Its
#' standard error uses the influence curve of the AUC in each fold (LeDell,
#' Petersen and van der Laan 2015). Pass `p_ref` to get the paired difference
#' `AUC(p) - AUC(p_ref)` on the same folds, with the standard error of the
#' difference of the influence curves.
#'
#' @param y Binary outcome.
#' @param p Held-out (cross-validated) predicted probabilities.
#' @param folds Fold label of every observation.
#' @param p_ref Optional reference predictions for a paired difference.
#' @param level Confidence level.
#' @return A list with `estimate`, `se`, `lower`, `upper` and the influence
#'   curve `ic` (one value per observation).
#' @export
cv_auc <- function(y, p, folds, p_ref = NULL, level = 0.95) {
  a <- .cv_auc_ic(y, p, folds)
  if (!is.null(p_ref)) {
    b <- .cv_auc_ic(y, p_ref, folds)
    a <- list(auc = a$auc - b$auc, ic = a$ic - b$ic)
  }
  se <- .ic_se(a$ic, folds)
  z <- stats::qnorm(1 - (1 - level) / 2)
  list(estimate = a$auc, se = se, lower = a$auc - z * se, upper = a$auc + z * se, ic = a$ic)
}

#' @noRd
.cv_auc_ic <- function(y, p, fold) {
  n <- length(y)
  w1 <- n / sum(y == 1); w0 <- n / sum(y == 0)
  ic <- numeric(n); fold_auc <- numeric(0)
  for (v in sort(unique(fold))) {
    i <- which(fold == v)
    cs <- y[i] == 1
    n1 <- sum(cs); n0 <- sum(!cs)
    r_all <- rank(p[i])
    share <- numeric(length(i))
    share[cs] <- (r_all[cs] - rank(p[i][cs])) / n0
    share[!cs] <- (n1 - (r_all[!cs] - rank(p[i][!cs]))) / n1
    auc_v <- mean(share[cs])
    ic[i] <- ifelse(cs, w1, w0) * (share - auc_v)
    fold_auc <- c(fold_auc, auc_v)
  }
  list(auc = mean(fold_auc), ic = ic)
}

#' @noRd
.ic_se <- function(ic, fold) sqrt(mean(tapply(ic^2, fold, mean)) / length(ic))
