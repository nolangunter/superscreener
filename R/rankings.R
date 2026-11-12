#' Ranking screeners
#'
#' Ranks every column of `X` by its importance for predicting the binary
#' outcome `Y` (1 = most important). These are the ranking screeners of
#' SuperScreener; the top `K` columns of a ranking form a candidate set.
#'
#' \describe{
#'   \item{`lrt`}{Univariate likelihood-ratio statistic: the drop in deviance
#'     when the column is added to an intercept-only logistic regression.}
#'   \item{`lasso`}{Order of entry along a logistic lasso path (200 penalty
#'     values, standardized columns). Columns that never enter share the last rank.}
#'   \item{`rf`}{Permutation importance of a `ranger` classification forest
#'     (`n_trees` trees).}
#'   \item{`consensus`}{Columns ordered by the unweighted mean of the `lrt`,
#'     `rf` and `lasso` ranks, ties broken by the `lrt` rank and then the `rf` rank.}
#'   \item{`spearman`}{Absolute Spearman correlation with the outcome, as
#'     `SuperLearner::screen.corRank(method = "spearman")`.}
#' }
#'
#' A column that is constant in `X` is placed last in the `lrt`, `lasso`, `rf`
#' and `consensus` rankings, and ties with other zero-correlation columns in
#' `spearman`. In the `lasso` ranking, columns that never enter the path share
#' the last rank, so [top_k()] fills a budget beyond the columns that entered
#' in column order.
#'
#' @param X Numeric data frame or matrix of predictors with column names.
#' @param Y Binary outcome coded 0/1.
#' @param rankings Which rankings to compute; any of [RANKINGS].
#' @param n_trees Number of trees for the random forest ranking.
#' @return A numeric matrix with one row per column of `X` (row names are the
#'   column names) and one column per ranking. Entry `[i, k]` is the rank of
#'   column `i` under ranking `k`.
#' @examples
#' set.seed(1)
#' X <- data.frame(matrix(rnorm(200 * 20), 200))
#' Y <- rbinom(200, 1, plogis(X$X1 - X$X2))
#' head(rank_columns(X, Y))
#' @export
rank_columns <- function(X, Y, rankings = RANKINGS, n_trees = 500) {
  X <- as.data.frame(X)
  rankings <- match.arg(rankings, RANKINGS, several.ok = TRUE)
  out <- matrix(NA_real_, ncol(X), length(rankings), dimnames = list(colnames(X), rankings))
  base <- intersect(rankings, c("lrt", "lasso", "rf", "consensus"))
  if (length(base)) out[, base] <- .rank_lrt_lasso_rf(X, Y, n_trees)[, base]
  if ("spearman" %in% rankings) out[, "spearman"] <- .rank_spearman(X, Y)
  out
}

#' Names of the ranking screeners
#' @export
RANKINGS <- c("lrt", "lasso", "rf", "consensus", "spearman")

#' The univariate LRT, random forest and lasso-entry ranks and their consensus.
#' Always computes all four, in this order (LRT, forest, lasso), so that the
#' random-number stream is the same whichever subset is requested.
#' @noRd
.rank_lrt_lasso_rf <- function(X, Y, n_trees = 500) {
  all_cols <- colnames(X)
  X <- X[, .nonconstant(X), drop = FALSE]
  p <- ncol(X)
  out <- matrix(length(all_cols), length(all_cols), 4,
                dimnames = list(all_cols, c("lrt", "lasso", "rf", "consensus")))
  if (p == 0) return(out)
  null_dev <- stats::glm.fit(matrix(1, length(Y), 1), Y, family = stats::binomial())$deviance
  lrt <- vapply(X, function(z) {
    fit <- tryCatch(suppressWarnings(stats::glm.fit(cbind(1, z), Y, family = stats::binomial())),
                    error = function(e) NULL)
    if (is.null(fit)) 0 else null_dev - fit$deviance
  }, numeric(1))
  rf <- ranger::ranger(x = X, y = factor(Y), num.trees = n_trees, importance = "permutation",
                       num.threads = 1, verbose = FALSE)
  las <- suppressWarnings(glmnet::glmnet(as.matrix(X), Y, family = "binomial", alpha = 1, standardize = TRUE, nlambda = 200))
  entry <- apply(as.matrix(las$beta != 0), 1, function(b) if (any(b)) which(b)[1] else NA_integer_)
  r_lrt <- rank(-lrt, ties.method = "min")
  r_rf <- rank(-rf$variable.importance[colnames(X)], ties.method = "min")
  r_lasso <- ifelse(is.na(entry), p, rank(entry, ties.method = "min", na.last = "keep"))
  ord <- order((r_lrt + r_rf + r_lasso) / 3, r_lrt, r_rf)
  r_cons <- integer(p); r_cons[ord] <- seq_len(p)
  out[colnames(X), ] <- cbind(r_lrt, r_lasso, r_rf, r_cons)
  out
}

#' @noRd
.rank_spearman <- function(X, Y) {
  r <- suppressWarnings(abs(stats::cor(as.matrix(X), Y, method = "spearman")[, 1]))
  r[!is.finite(r)] <- 0
  stats::setNames(rank(-r, ties.method = "first"), colnames(X))
}

#' Top K columns of a ranking
#'
#' @param ranks Named numeric vector of ranks (1 = most important), e.g. one
#'   column of [rank_columns()].
#' @param K Number of columns to keep.
#' @return Character vector of the `min(K, length(ranks))` best-ranked column names.
#' @export
top_k <- function(ranks, K) {
  names(sort(ranks, decreasing = FALSE))[seq_len(min(K, length(ranks)))]
}

#' @noRd
.nonconstant <- function(X) vapply(as.data.frame(X), function(z) length(unique(z)) > 1, logical(1))
