#' Fold assignment stratified by outcome
#'
#' Splits the cases and the non-cases separately into `V` folds of nearly equal
#' size, so every fold has about the same proportion of cases.
#' @param Y Binary outcome.
#' @param V Number of folds.
#' @return Integer vector of fold labels `1..V`.
#' @export
stratified_folds <- function(Y, V = 5) {
  f <- integer(length(Y))
  for (cl in unique(Y)) {
    idx <- which(Y == cl)
    f[idx] <- sample(rep_len(seq_len(V), length(idx)))
  }
  f
}

#' Expand categorical variables into indicator columns
#'
#' SuperScreener charges its price per *variable*, and a categorical variable
#' coded as several indicator columns costs one unit. This helper turns a data
#' frame with factor, character or logical columns into a numeric data frame
#' and returns the map from each column to its variable, ready for the
#' `groups` argument of [superscreener()].
#'
#' Numeric and logical columns are kept as they are (logical as 0/1). A factor
#' or character variable with `L` levels becomes `L - 1` indicator columns
#' named `variable_level` (first level is the reference). Missing values are
#' not allowed: impute them (and add missingness indicators if you like) first.
#'
#' @param data A data frame of predictors (without the outcome).
#' @return A list with `X` (numeric data frame) and `groups` (named character
#'   vector: names are the columns of `X`, values the original variables).
#' @examples
#' d <- data.frame(age = c(30, 45, 60, 52), sex = c("f", "m", "f", "m"),
#'                 injury = factor(c("head", "limb", "torso", "limb")))
#' expand_factors(d)
#' @export
expand_factors <- function(data) {
  data <- as.data.frame(data)
  if (anyNA(data)) stop("`data` has missing values; impute them before calling expand_factors().")
  cols <- list(); groups <- character(0)
  for (v in names(data)) {
    x <- data[[v]]
    if (is.logical(x)) x <- as.numeric(x)
    if (is.character(x)) x <- factor(x)
    if (is.factor(x)) {
      x <- droplevels(x)
      if (nlevels(x) < 2) { cols[[v]] <- rep(0, length(x)); groups[v] <- v; next }
      mm <- stats::model.matrix(~ x)[, -1, drop = FALSE]
      nm <- make.names(paste0(v, "_", levels(x)[-1]), unique = TRUE)
      for (k in seq_len(ncol(mm))) { cols[[nm[k]]] <- unname(mm[, k]); groups[nm[k]] <- v }
    } else if (is.numeric(x)) {
      cols[[v]] <- x; groups[v] <- v
    } else stop(sprintf("Column `%s` has unsupported type %s.", v, class(x)[1]))
  }
  all_names <- unlist(lapply(names(data), function(v) {
    x <- data[[v]]
    if (is.character(x) || is.factor(x)) {
      lv <- levels(droplevels(factor(x)))
      if (length(lv) < 2) v else make.names(paste0(v, "_", lv[-1]), unique = TRUE)
    } else v
  }))
  dup <- unique(all_names[duplicated(all_names)])
  if (length(dup)) stop("Indicator column name(s) collide with existing columns: ", paste(dup, collapse = ", "),
                        ". Rename those columns first.")
  list(X = as.data.frame(cols, check.names = FALSE), groups = groups)
}

#' Check and standardize X and Y. Returns list(X = data frame, Y = numeric 0/1).
#' `V` is the number of folds that will be made: each class needs at least V observations.
#' @noRd
.check_xy <- function(X, Y, V = 2) {
  if (!is.data.frame(X) && !is.matrix(X)) stop("`X` must be a data frame or matrix.")
  X <- as.data.frame(X)
  if (is.null(colnames(X)) || anyDuplicated(colnames(X))) stop("`X` needs unique column names.")
  bad <- !vapply(X, is.numeric, logical(1))
  if (any(bad)) stop("Non-numeric columns in `X` (", paste(utils::head(names(X)[bad], 5), collapse = ", "),
                     "). Use expand_factors() first.")
  if (anyNA(X)) stop("`X` has missing values; impute them first.")
  if (length(Y) != nrow(X)) stop("`Y` and `X` have different numbers of rows.")
  if (is.factor(Y)) Y <- as.character(Y)
  if (is.logical(Y)) Y <- as.numeric(Y)
  Y <- suppressWarnings(as.numeric(Y))
  if (anyNA(Y) || !all(Y %in% c(0, 1))) stop("`Y` must be coded 0/1 (or TRUE/FALSE) with no missing values.")
  if (length(unique(Y)) < 2) stop("`Y` has only one class.")
  if (min(table(Y)) < V)
    stop(sprintf("The rarer outcome class has %d observations, fewer than the %d folds requested.",
                 min(table(Y)), V))
  list(X = X, Y = Y)
}

#' @noRd
.check_args <- function(K, rankings, self_sized, price, p) {
  bad_r <- setdiff(rankings, RANKINGS); bad_s <- setdiff(self_sized, SELF_SIZED)
  if (length(bad_r)) stop("Unknown ranking(s): ", paste(bad_r, collapse = ", "), ". Choose from ",
                          paste(RANKINGS, collapse = ", "), ".")
  if (length(bad_s)) stop("Unknown self-sized screener(s): ", paste(bad_s, collapse = ", "), ". Choose from ",
                          paste(SELF_SIZED, collapse = ", "), ".")
  if (length(rankings) && (!length(K) || any(!is.finite(K)) || any(K < 1) || any(K != round(K))))
    stop("`K` must be positive whole numbers.")
  if (any(!is.finite(price)) || any(price < 0)) stop("Prices must be >= 0.")
  invisible(TRUE)
}
