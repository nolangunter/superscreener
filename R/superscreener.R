#' SuperScreener: choose a screener and the number of variables by
#' cross-validated risk plus a price per variable
#'
#' Runs the six steps of the SuperScreener algorithm on a binary outcome:
#' \enumerate{
#'   \item Split the sample into `V` folds stratified by outcome.
#'   \item In each training sample, run every screener, take the top `K`
#'     columns of each ranking for every `K` and the set of each self-sized
#'     screener, and fit the same Super Learner to each set.
#'   \item Predict the held-out fold with every candidate. The cross-validated
#'     risk \eqn{\hat R(j)} of candidate \eqn{j} is the mean loss of its held-out
#'     predictions; \eqn{\hat R(\emptyset)} is that of the intercept-only model.
#'   \item Run every screener on the whole sample; \eqn{s_j} is the cost (by
#'     default the number of variables) of the set candidate \eqn{j} selects there.
#'   \item Set \eqn{\lambda = c\{\hat R(\emptyset) - \min_j \hat R(j)\}} (0 if no
#'     candidate beats the intercept) and choose the candidate minimizing
#'     \eqn{\hat R(j) + \lambda s_j}, ties going to the smaller set.
#'   \item Fit the Super Learner to the chosen candidate's whole-sample set.
#' }
#' The relative price `price` (\eqn{c}) is the share of the best candidate's
#' risk reduction over the intercept-only model that one more variable has to
#' buy. `price = 0` is the discrete Super Learner over the candidates.
#'
#' The cross-validated AUC stored in `candidates` is optimistic for the chosen
#' candidate, because the same folds chose it. Use [nested_cv_superscreener()]
#' for an honest estimate of the accuracy of the whole procedure.
#'
#' @param X Numeric data frame or matrix of predictors with unique column
#'   names and no missing values. Use [expand_factors()] for categorical data.
#' @param Y Binary outcome coded 0/1.
#' @param K Budgets (numbers of columns) at which each ranking is used.
#' @param rankings Ranking screeners to include; see [rank_columns()].
#'   `character(0)` for none.
#' @param self_sized Self-sized screeners to include; see [self_sized_sets()].
#'   `character(0)` for none.
#' @param price Relative price per variable, \eqn{c}. The chapter uses 0.01.
#' @param V Number of cross-validation folds for the selector.
#' @param loss Selection loss: `"brier"` (squared error, the default) or
#'   `"neg_loglik"`.
#' @param sl_library Super Learner library for the final model.
#' @param selection_library Super Learner library fit to every candidate inside
#'   the cross-validation. Defaults to `sl_library`; a smaller library makes
#'   the selector faster.
#' @param sl_folds Folds of the Super Learner's own internal cross-validation.
#' @param groups Optional named character vector mapping each column of `X`
#'   (names) to the variable it belongs to (values), e.g. from
#'   [expand_factors()]. A variable is charged once however many of its
#'   columns are selected. Default: every column is its own variable.
#' @param cost Optional named numeric vector of the cost of each variable
#'   (names = variables). Default 1 for every variable, so \eqn{s_j} is the
#'   number of variables.
#' @param min_size Drop candidates whose cost is below this before choosing
#'   (the chapter's post hoc "at least five variables" rule is `min_size = 5`).
#' @param n_trees Trees for the random forest ranking.
#' @param seed Optional random seed, for reproducible folds and fits.
#' @param keep_data Keep `X` and `Y` in the result, so that [set_price()]
#'   can refit at another price without rerunning the cross-validation.
#' @param verbose Print progress.
#' @return An object of class `"superscreener"`, a list with
#'   \item{selected}{columns used by the final model}
#'   \item{variables}{variables (groups) used by the final model}
#'   \item{chosen}{the chosen row of `candidates`}
#'   \item{candidates}{one row per candidate: `candidate`, `screener`, `K`
#'     (`NA` for self-sized), `n_columns`, `size` (\eqn{s_j}), `cv_risk`,
#'     `cv_auc`, `cv_auc_se`, `penalized_risk`}
#'   \item{sets}{whole-sample set of every candidate}
#'   \item{risk_null, lambda, price, loss}{selection quantities}
#'   \item{cv_predictions, folds}{held-out predictions of every candidate}
#'   \item{fit}{the final Super Learner (see [predict.superscreener()])}
#' @examples
#' \donttest{
#' design <- make_design(p = 60)
#' dat <- simulate_data(300, design, seed = 1)
#' fit <- superscreener(dat$X, dat$Y, K = c(5, 10, 20), price = 0.01,
#'                      sl_library = c("SL.glm", "SL.mean"), seed = 1)
#' fit
#' price_path(fit)
#' }
#' @export
superscreener <- function(X, Y, K = seq(5, 30, by = 5), rankings = RANKINGS, self_sized = SELF_SIZED,
                          price = 0.01, V = 5, loss = c("brier", "neg_loglik"),
                          sl_library = sl_library_default(), selection_library = sl_library, sl_folds = 5,
                          groups = NULL, cost = NULL, min_size = 0, n_trees = 500,
                          seed = NULL, keep_data = TRUE, verbose = TRUE) {
  loss <- match.arg(loss)
  chk <- .check_xy(X, Y, V); X <- chk$X; Y <- chk$Y
  .check_args(K, rankings, self_sized, price, ncol(X))
  if (!is.null(seed)) set.seed(seed)
  cv <- .cv_candidates(X, Y, K = K, rankings = rankings, self_sized = self_sized, V = V, loss = loss,
                       library = selection_library, sl_folds = sl_folds, groups = groups, cost = cost,
                       n_trees = n_trees, verbose = verbose)
  out <- c(cv, list(price = price, min_size = min_size, sl_library = sl_library, sl_folds = sl_folds))
  if (keep_data) { out$X <- X; out$Y <- Y }
  class(out) <- "superscreener"
  .finalize(out, X, Y, price, min_size, verbose)
}

#' Steps 1 to 4: candidates, their cross-validated risks and whole-sample sets.
#' @noRd
.cv_candidates <- function(X, Y, K, rankings, self_sized, V, loss, library, sl_folds, groups, cost,
                           n_trees, verbose = FALSE) {
  cols <- colnames(X)
  groups <- .resolve_groups(groups, cols)
  cost <- .resolve_cost(cost, groups)
  ## a budget above the number of columns would repeat the full set under another label
  K <- sort(unique(pmin(as.integer(K), length(cols))))
  cand <- rbind(
    if (length(rankings)) data.frame(screener = rep(rankings, each = length(K)), K = rep(K, length(rankings))),
    if (length(self_sized)) data.frame(screener = self_sized, K = NA_integer_))
  if (is.null(cand) || !nrow(cand)) stop("No candidates: give at least one ranking or self-sized screener.")
  cand$candidate <- ifelse(is.na(cand$K), cand$screener, paste0(cand$screener, "_K", cand$K))
  sets_of <- function(Xs, Ys) {
    rk <- if (length(rankings)) rank_columns(Xs, Ys, rankings, n_trees)
    sz <- if (length(self_sized)) self_sized_sets(Xs, Ys, self_sized)
    stats::setNames(lapply(seq_len(nrow(cand)), function(j)
      if (is.na(cand$K[j])) sz[[cand$screener[j]]] else top_k(rk[, cand$screener[j]], cand$K[j])),
      cand$candidate)
  }
  lossfun <- if (loss == "brier") brier else neg_loglik

  ## steps 1-3
  folds <- stratified_folds(Y, V)
  pred <- matrix(NA_real_, length(Y), nrow(cand), dimnames = list(NULL, cand$candidate))
  pred_null <- numeric(length(Y))
  for (v in seq_len(V)) {
    if (verbose) message(sprintf("fold %d of %d: screening and fitting %d candidates", v, V, nrow(cand)))
    trn <- folds != v; tst <- !trn
    s_v <- sets_of(X[trn, , drop = FALSE], Y[trn])
    pred_null[tst] <- mean(Y[trn])
    for (j in seq_len(nrow(cand))) {
      cc <- s_v[[j]]
      pred[tst, j] <- fit_sl(X[trn, cc, drop = FALSE], Y[trn], X[tst, cc, drop = FALSE], library, sl_folds)
    }
  }
  cand$cv_risk <- apply(pred, 2, function(p) lossfun(Y, p))
  au <- lapply(seq_len(nrow(cand)), function(j) cv_auc(Y, pred[, j], folds))
  cand$cv_auc <- vapply(au, `[[`, 0, "estimate")
  cand$cv_auc_se <- vapply(au, `[[`, 0, "se")

  ## step 4
  if (verbose) message("screening the whole sample")
  sets <- sets_of(X, Y)
  cand$n_columns <- lengths(sets)
  cand$size <- vapply(sets, function(s) sum(cost[unique(groups[s])]), numeric(1))
  list(candidates = cand[, c("candidate", "screener", "K", "n_columns", "size", "cv_risk", "cv_auc", "cv_auc_se")],
       sets = sets, risk_null = lossfun(Y, pred_null), loss = loss, groups = groups, cost = cost,
       cv_predictions = pred, folds = folds)
}

#' Step 5: choose the candidate with the smallest penalized risk
#'
#' @param candidates Data frame with columns `cv_risk`, `size` and `K`
#'   (`NA` for self-sized screeners), e.g. `fit$candidates`.
#' @param price Relative price \eqn{c}.
#' @param risk_null Cross-validated risk of the intercept-only model.
#' @param min_size Candidates with `size < min_size` are dropped first.
#' @return The chosen row of `candidates`, with `lambda` and `penalized_risk`
#'   added. The price is computed from the candidates left after `min_size`.
#' @export
choose_candidate <- function(candidates, price, risk_null, min_size = 0) {
  cand <- candidates[candidates$size >= min_size, , drop = FALSE]
  if (!nrow(cand)) stop("No candidate has size >= min_size.")
  lam <- price * max(risk_null - min(cand$cv_risk), 0)
  cand$penalized_risk <- cand$cv_risk + lam * cand$size
  k_tie <- ifelse(is.na(cand$K), 0, cand$K)
  out <- cand[order(cand$penalized_risk, cand$size, k_tie)[1], , drop = FALSE]
  out$lambda <- lam
  out
}

#' Choices across a range of prices
#'
#' Applies step 5 at each price without refitting anything. Useful for showing
#' how many variables are given up for how much accuracy.
#' @param object A `"superscreener"` fit.
#' @param prices Relative prices \eqn{c}.
#' @param min_size As in [superscreener()].
#' @return A data frame with one row per price.
#' @export
price_path <- function(object, prices = c(0, 0.0025, 0.005, 0.01, 0.02, 0.04), min_size = object$min_size) {
  do.call(rbind, lapply(prices, function(c) {
    ch <- choose_candidate(object$candidates, c, object$risk_null, min_size)
    data.frame(price = c, lambda = ch$lambda, candidate = ch$candidate, size = ch$size,
               n_columns = ch$n_columns, cv_risk = ch$cv_risk, cv_auc = ch$cv_auc, row.names = NULL)
  }))
}

#' Refit a SuperScreener at another price
#'
#' Re-chooses the candidate at `price` from the stored cross-validated risks
#' and refits the final Super Learner. Needs `keep_data = TRUE`.
#' @param object A `"superscreener"` fit.
#' @param price New relative price.
#' @param min_size As in [superscreener()].
#' @param verbose Print progress.
#' @return An updated `"superscreener"` object.
#' @export
set_price <- function(object, price, min_size = object$min_size, verbose = FALSE) {
  if (is.null(object$X)) stop("Refit needs the data: rerun superscreener() with keep_data = TRUE.")
  .finalize(object, object$X, object$Y, price, min_size, verbose)
}

#' Steps 5 and 6.
#' @noRd
.finalize <- function(object, X, Y, price, min_size, verbose) {
  ch <- choose_candidate(object$candidates, price, object$risk_null, min_size)
  cand <- object$candidates
  lam <- ch$lambda
  cand$penalized_risk <- cand$cv_risk + lam * cand$size
  object$candidates <- cand
  object$chosen <- ch
  object$price <- price; object$min_size <- min_size; object$lambda <- lam
  object$selected <- object$sets[[ch$candidate]]
  object$variables <- unique(unname(object$groups[object$selected]))
  if (verbose) message(sprintf("fitting the final Super Learner on %d columns", length(object$selected)))
  object$fit <- fit_sl(X[, object$selected, drop = FALSE], Y, X[1, object$selected, drop = FALSE],
                       object$sl_library, object$sl_folds, return_fit = TRUE)
  object$fit$pred <- NULL
  failed <- if (!is.null(object$fit$sl)) object$fit$sl$libraryNames[object$fit$sl$errorsInLibrary] else character(0)
  if (length(failed))
    warning("These learners failed in the final Super Learner and were given weight 0: ",
            paste(failed, collapse = ", "), call. = FALSE)
  object
}

#' Predict from a SuperScreener fit
#' @param object A `"superscreener"` fit.
#' @param newdata Data frame containing (at least) the selected columns.
#' @param ... Unused.
#' @return Predicted probabilities.
#' @export
predict.superscreener <- function(object, newdata, ...) {
  newdata <- as.data.frame(newdata)
  if (is.null(object$fit$sl)) return(rep(object$fit$prior, nrow(newdata)))
  miss <- setdiff(object$fit$columns, colnames(newdata))
  if (length(miss)) stop("`newdata` is missing columns: ", paste(utils::head(miss, 5), collapse = ", "))
  as.numeric(stats::predict(object$fit$sl, newdata = newdata[, object$fit$columns, drop = FALSE],
                            onlySL = TRUE)$pred)
}

#' @export
print.superscreener <- function(x, ...) {
  ch <- x$chosen
  cat("SuperScreener\n")
  cat(sprintf("  candidates: %d (%s loss, %d-fold CV)\n", nrow(x$candidates), x$loss, length(unique(x$folds))))
  cat(sprintf("  price c = %g  ->  lambda = %.3g per unit of cost\n", x$price, x$lambda))
  cat(sprintf("  chosen: %s  (%d columns, %d variables, cost %g)\n", ch$candidate, length(x$selected),
              length(x$variables), ch$size))
  cat(sprintf("  CV risk %.4f (intercept only %.4f); CV AUC %.3f (optimistic: chosen on these folds)\n",
              ch$cv_risk, x$risk_null, ch$cv_auc))
  shown <- utils::head(x$variables, 15)
  cat("  variables:", paste(shown, collapse = ", "), if (length(x$variables) > 15) "...", "\n")
  invisible(x)
}

#' @noRd
.resolve_groups <- function(groups, cols) {
  if (is.null(groups)) return(stats::setNames(cols, cols))
  miss <- setdiff(cols, names(groups))
  if (length(miss)) stop("`groups` has no entry for columns: ", paste(utils::head(miss, 5), collapse = ", "))
  ## as.character: a factor would index `cost` by its integer codes instead of its labels
  stats::setNames(as.character(groups[cols]), cols)
}

#' @noRd
.resolve_cost <- function(cost, groups) {
  vars <- unique(unname(groups))
  if (is.null(cost)) return(stats::setNames(rep(1, length(vars)), vars))
  miss <- setdiff(vars, names(cost))
  if (length(miss)) stop("`cost` has no entry for variables: ", paste(utils::head(miss, 5), collapse = ", "))
  if (any(!is.finite(cost[vars])) || any(cost[vars] < 0)) stop("`cost` must be finite and >= 0.")
  cost[vars]
}
