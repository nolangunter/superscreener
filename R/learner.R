#' The default Super Learner library
#'
#' A random forest with 300 trees ([SL.ranger_300]), a logistic lasso
#' (`SL.glmnet`), gradient boosting with 200 trees ([SL.xgboost_200]) and the
#' proportion of cases (`SL.mean`). This is the downstream learner used in the
#' simulation study. Any `SuperLearner` library can be used instead.
#' @return A character vector of `SuperLearner` wrapper names.
#' @export
sl_library_default <- function() c("SL.ranger_300", "SL.glmnet", "SL.xgboost_200", "SL.mean")

#' Learner wrappers with fixed tree counts
#'
#' `SuperLearner::SL.ranger` with 300 trees and `SuperLearner::SL.xgboost`
#' with 200 trees.
#' @param ... Passed to the `SuperLearner` wrapper.
#' @export
SL.ranger_300 <- function(...) SuperLearner::SL.ranger(..., num.trees = 300)

#' @rdname SL.ranger_300
#' @export
SL.xgboost_200 <- function(...) SuperLearner::SL.xgboost(..., ntrees = 200)

#' Environment in which SuperLearner looks up learner and screener names: this
#' package's wrappers, then SuperLearner's own (SL.glmnet, All, ...), then the
#' global environment and search path (user-defined learners).
#' @noRd
.sl_env <- function() {
  e <- new.env(parent = asNamespace("SuperLearner"))
  e$SL.ranger_300 <- SL.ranger_300
  e$SL.xgboost_200 <- SL.xgboost_200
  e
}

#' Fit a Super Learner to a set of columns and predict new data
#'
#' A binomial Super Learner whose weights minimize the cross-validated negative
#' log-likelihood (`method.NNloglik`), with folds stratified by outcome.
#' Columns that are constant in the training data are dropped. A learner that
#' fails on the set (for example `SL.glmnet` with a single column) gets weight 0,
#' as in `SuperLearner`. If no column is left, or the whole fit fails, the
#' prediction is the proportion of cases in the training data.
#'
#' @param Xtr,Ytr Training predictors (data frame) and binary outcome.
#' @param Xte Data frame to predict, with the same columns as `Xtr`.
#' @param library A `SuperLearner` library.
#' @param V Number of folds for the Super Learner's internal cross-validation.
#' @param return_fit If `TRUE`, return a list with the predictions and the
#'   fitted object (used for the final model) instead of the predictions alone.
#' @return Numeric vector of predicted probabilities for the rows of `Xte`
#'   (or a list, see `return_fit`).
#' @export
fit_sl <- function(Xtr, Ytr, Xte, library = sl_library_default(), V = 5, return_fit = FALSE) {
  Xtr <- as.data.frame(Xtr); Xte <- as.data.frame(Xte)
  if (ncol(Xtr) > 0) {
    Xtr <- Xtr[, .nonconstant(Xtr), drop = FALSE]
    Xte <- Xte[, colnames(Xtr), drop = FALSE]
  }
  prior <- mean(Ytr)
  fallback <- function() {
    pr <- rep(prior, nrow(Xte))
    if (return_fit) list(pred = pr, sl = NULL, columns = character(0), prior = prior) else pr
  }
  if (ncol(Xtr) == 0) return(fallback())
  ## A learner that fails on this set (e.g. SL.glmnet with a single column) is dropped from the
  ## ensemble by SuperLearner, which also prints the learner's error; keep the console quiet.
  old <- options(try.outFile = nullfile()); on.exit(options(old), add = TRUE)
  sl <- tryCatch(suppressWarnings(SuperLearner::SuperLearner(
    Y = Ytr, X = Xtr, newX = Xte, family = stats::binomial(), SL.library = library,
    method = "method.NNloglik", cvControl = list(V = V, stratifyCV = TRUE), env = .sl_env())),
    error = function(e) { warning("Super Learner fit failed, predicting the case proportion: ",
                                  conditionMessage(e), call. = FALSE); NULL })
  ## If every learner that could fit was SL.mean, method.NNloglik gives all weights 0 and
  ## SuperLearner predicts 0 for everyone; use the case proportion instead.
  if (is.null(sl) || !any(sl$coef > 0)) return(fallback())
  pr <- as.numeric(sl$SL.predict)
  if (return_fit) list(pred = pr, sl = sl, columns = colnames(Xtr), prior = prior) else pr
}
