#' Self-sized screeners
#'
#' Screeners that choose their own number of columns. Each returns a set of
#' column names, which may be empty (`lasso_1se`, `stabsel`).
#'
#' \describe{
#'   \item{`sl_glmnet`}{`SuperLearner::screen.glmnet`: nonzero coefficients of a
#'     logistic lasso at the penalty minimizing 10-fold cross-validated deviance;
#'     keeps at least two columns.}
#'   \item{`sl_corP`}{`SuperLearner::screen.corP`: Pearson correlation test
#'     p-value at most 0.1; keeps at least the two smallest p-values.}
#'   \item{`lasso_1se`}{Nonzero logistic lasso coefficients at `lambda.1se`
#'     (10-fold CV). No minimum, so it can keep no columns.}
#'   \item{`stabsel`}{Stability selection with complementary pairs (Shah and
#'     Samworth 2013): in each of `2 * stabsel_pairs` half-samples, take the first
#'     `q = max(2, round(sqrt(0.5 p)))` columns to enter the lasso path; keep the columns
#'     selected in at least 75\% of half-samples. With this `q` the
#'     Meinshausen-Buhlmann bound on expected false selections is about 1.}
#' }
#'
#' @inheritParams rank_columns
#' @param screeners Which screeners to run; any of [SELF_SIZED].
#' @param stabsel_pairs Number of complementary pairs for stability selection.
#' @param seed Optional seed set before the screeners run.
#' @return A named list of character vectors, one per screener.
#' @export
self_sized_sets <- function(X, Y, screeners = SELF_SIZED, stabsel_pairs = 25, seed = NULL) {
  screeners <- match.arg(screeners, SELF_SIZED, several.ok = TRUE)
  if (!is.null(seed)) set.seed(seed)
  X <- as.data.frame(X)
  Xk <- X[, .nonconstant(X), drop = FALSE]
  Xm <- as.matrix(Xk)
  ## the lasso-based screeners need at least two non-constant columns; with fewer, keep them all
  if (ncol(Xm) < 2) return(stats::setNames(rep(list(colnames(Xk)), length(screeners)), screeners))
  out <- list()
  if ("sl_glmnet" %in% screeners)
    out$sl_glmnet <- colnames(Xk)[SuperLearner::screen.glmnet(Y, Xk, family = stats::binomial())]
  if ("sl_corP" %in% screeners)
    out$sl_corP <- colnames(Xk)[suppressWarnings(SuperLearner::screen.corP(Y, Xk, family = stats::binomial()))]
  if ("lasso_1se" %in% screeners) {
    cvf <- glmnet::cv.glmnet(Xm, Y, family = "binomial", nfolds = 10)
    b <- as.matrix(stats::coef(cvf, s = "lambda.1se"))[-1, 1]
    out$lasso_1se <- names(b)[b != 0]
  }
  if ("stabsel" %in% screeners) out$stabsel <- .stability_selection(Xm, Y, stabsel_pairs)
  out[screeners]
}

#' Names of the self-sized screeners
#' @export
SELF_SIZED <- c("sl_glmnet", "sl_corP", "lasso_1se", "stabsel")

#' @noRd
.stability_selection <- function(Xm, Y, pairs = 25, threshold = 0.75) {
  n <- nrow(Xm); p <- ncol(Xm); q <- max(2, round(sqrt(0.5 * p)))
  hits <- stats::setNames(numeric(p), colnames(Xm))
  for (b in seq_len(pairs)) {
    half <- sample(n, floor(n / 2))
    for (idx in list(half, setdiff(seq_len(n), half))) {
      # a half-sample with fewer than 2 observations in a class cannot be fit; it selects nothing
      if (min(table(factor(Y[idx], levels = 0:1))) < 2) next
      fit <- glmnet::glmnet(Xm[idx, , drop = FALSE], Y[idx], family = "binomial", nlambda = 100, dfmax = q + 5)
      entry <- apply(as.matrix(fit$beta != 0), 1, function(z) if (any(z)) which(z)[1] else Inf)
      first <- names(sort(entry))[seq_len(min(q, sum(is.finite(entry))))]
      hits[first] <- hits[first] + 1
    }
  }
  names(hits)[hits / (2 * pairs) >= threshold]
}
