#' Nested cross-validation of SuperScreener
#'
#' Reports the accuracy and stability of the whole SuperScreener procedure,
#' including its choice of screener and number of variables. In each training
#' sample of an outer `V_outer`-fold cross-validation, SuperScreener runs all
#' its steps (inner cross-validation, choice at each price, final fit), and the
#' resulting model predicts the held-out fold. The held-out predictions give
#' the nested cross-validated AUC, with an influence-curve interval (LeDell et
#' al. 2015), and the paired difference from the model with every column.
#'
#' The inner cross-validation is run once per outer training sample and reused
#' for every price, so a grid of prices costs little more than one price.
#' With `repeats > 1`, estimates and standard errors are averaged over the
#' repeated splits (the interval does not account for this averaging).
#'
#' @inheritParams superscreener
#' @param prices Relative prices \eqn{c} at which to evaluate the procedure.
#' @param V_outer Number of outer folds.
#' @param repeats Number of repeated outer splits.
#' @param include_full Also fit the model with every column, for reference.
#' @param cores Outer folds to run in parallel, in separate R processes (needs the
#'   package installed). Results do not depend on `cores`.
#' @param seed Seed. Each outer fold uses its own seed derived from it.
#' @param ... Further arguments to [superscreener()] (e.g. `K`, `rankings`,
#'   `self_sized`, `V`, `loss`, `selection_library`, `min_size`).
#' @return An object of class `"nested_superscreener"`, a list with
#'   \item{summary}{one row per price (and `"full"`): mean number of
#'     variables and columns, nested cvAUC with 95\% interval, difference
#'     from the full model with interval, stability index of the selected
#'     variables across the outer training samples}
#'   \item{frequency}{selection frequency of every variable at every price}
#'   \item{choices}{the candidate chosen in every outer fold at every price}
#'   \item{predictions}{held-out predictions, one column per price}
#' @examples
#' \donttest{
#' design <- make_design(p = 60)
#' dat <- simulate_data(300, design, seed = 1)
#' ncv <- nested_cv_superscreener(dat$X, dat$Y, prices = c(0, 0.01, 0.04),
#'                                K = c(5, 10, 20), sl_library = c("SL.glm", "SL.mean"),
#'                                V_outer = 3, seed = 1)
#' ncv$summary
#' }
#' @export
nested_cv_superscreener <- function(X, Y, prices = c(0, 0.0025, 0.005, 0.01, 0.02, 0.04),
                                    V_outer = 5, repeats = 1, include_full = TRUE,
                                    sl_library = sl_library_default(), sl_folds = 5,
                                    groups = NULL, cost = NULL, cores = 1, seed = 1,
                                    verbose = TRUE, ...) {
  ## evaluate every argument now: parallel workers receive one_fold's environment and cannot
  ## evaluate promises that refer to the caller's variables
  for (a in c("prices", "V_outer", "repeats", "include_full", "sl_library", "sl_folds", "cost", "seed", "verbose"))
    force(get(a))
  dots <- list(...)
  inner_V <- if (is.null(dots$V)) 5 else dots$V
  ## each outer training sample must still have enough cases for the inner folds
  chk <- .check_xy(X, Y, max(V_outer, ceiling(inner_V * V_outer / (V_outer - 1)))); X <- chk$X; Y <- chk$Y
  .check_args(if (is.null(dots$K)) 5 else dots$K, if (is.null(dots$rankings)) RANKINGS else dots$rankings,
              if (is.null(dots$self_sized)) SELF_SIZED else dots$self_sized, prices, ncol(X))
  if (is.null(seed)) seed <- sample.int(1e6, 1)
  min_size <- if (is.null(dots$min_size)) 0 else dots$min_size
  dots$min_size <- NULL
  groups <- .resolve_groups(groups, colnames(X))
  vars <- unique(unname(groups))
  set.seed(seed)
  outer <- replicate(repeats, stratified_folds(Y, V_outer))
  tasks <- expand.grid(v = seq_len(V_outer), r = seq_len(repeats))

  one_fold <- function(t) {
    r <- tasks$r[t]; v <- tasks$v[t]
    set.seed(seed + 1000 * r + v)
    trn <- outer[, r] != v; tst <- !trn
    if (verbose) message(sprintf("repeat %d, outer fold %d", r, v))
    inner <- do.call(.cv_candidates, c(list(X = X[trn, , drop = FALSE], Y = Y[trn], groups = groups,
                                            cost = cost, verbose = FALSE),
                                       .inner_defaults(dots, sl_library, sl_folds)))
    cache <- list(); res <- list()
    for (c in prices) {
      ch <- choose_candidate(inner$candidates, c, inner$risk_null, min_size)
      cols <- inner$sets[[ch$candidate]]
      key <- paste0("set:", paste(sort(cols), collapse = "|"))   # never "", which [[ ]] cannot store
      if (is.null(cache[[key]]))
        cache[[key]] <- fit_sl(X[trn, cols, drop = FALSE], Y[trn], X[tst, cols, drop = FALSE], sl_library, sl_folds)
      res[[length(res) + 1]] <- list(price = c, candidate = ch$candidate, columns = cols,
                                     variables = unique(unname(groups[cols])), pred = cache[[key]])
    }
    full <- if (include_full) fit_sl(X[trn, , drop = FALSE], Y[trn], X[tst, , drop = FALSE], sl_library, sl_folds)
    list(r = r, v = v, tst = which(tst), res = res, full = full)
  }
  ## Parallel folds run in fresh R processes (a socket cluster), not forks: on macOS, forked
  ## workers that run ranger and then glmnet can slow down by orders of magnitude.
  runs <- if (cores > 1) {
    cl <- parallel::makeCluster(min(cores, nrow(tasks)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
      suppressPackageStartupMessages(library(superscreener))
    })
    ## learners the user defined in their session (e.g. SL.mine) do not exist in the workers
    sel_lib <- if (is.null(dots$selection_library)) sl_library else dots$selection_library
    learners <- unique(unlist(c(sl_library, sel_lib)))
    own <- learners[vapply(learners, exists, logical(1), envir = globalenv(), inherits = FALSE)]
    if (length(own)) parallel::clusterExport(cl, own, envir = globalenv())
    not_found <- function(l) l[!vapply(l, function(f)
      exists(f, envir = new.env(parent = asNamespace("SuperLearner"))) ||
        exists(f, envir = asNamespace("superscreener")), logical(1))]
    environment(not_found) <- globalenv()      # do not ship this function's frame (with X) to the worker
    missing <- unlist(parallel::clusterCall(cl[1], not_found, learners))
    if (length(missing)) stop("Learner(s) not found in the parallel workers: ", paste(missing, collapse = ", "),
                              ". Define them in the global environment or use cores = 1.")
    parallel::parLapplyLB(cl, seq_len(nrow(tasks)), one_fold)
  } else lapply(seq_len(nrow(tasks)), one_fold)

  ## held-out predictions: [observation, price, repeat]
  P <- array(NA_real_, c(length(Y), length(prices), repeats))
  Pfull <- matrix(NA_real_, length(Y), repeats)
  choices <- list()
  for (run in runs) {
    for (k in seq_along(prices)) {
      P[run$tst, k, run$r] <- run$res[[k]]$pred
      choices[[length(choices) + 1]] <- data.frame(repeat_id = run$r, fold = run$v, price = prices[k],
                                                   candidate = run$res[[k]]$candidate,
                                                   n_columns = length(run$res[[k]]$columns),
                                                   n_variables = length(run$res[[k]]$variables))
    }
    if (include_full) Pfull[run$tst, run$r] <- run$full
  }
  choices <- do.call(rbind, choices)
  M <- length(runs)
  freq <- sapply(seq_along(prices), function(k)
    selection_frequency(lapply(runs, function(run) run$res[[k]]$variables), vars))
  freq <- matrix(freq, nrow = length(vars), dimnames = list(vars, paste0("c=", prices)))

  avg <- function(f) {           # average estimate and SE over repeats
    z <- lapply(seq_len(repeats), function(r) f(r))
    c(estimate = mean(vapply(z, `[[`, 0, "estimate")), se = mean(vapply(z, `[[`, 0, "se")))
  }
  rows <- lapply(seq_along(prices), function(k) {
    a <- avg(function(r) cv_auc(Y, P[, k, r], outer[, r]))
    d <- if (include_full) avg(function(r) cv_auc(Y, P[, k, r], outer[, r], p_ref = Pfull[, r])) else c(estimate = NA, se = NA)
    ch <- choices[choices$price == prices[k], ]
    data.frame(method = "superscreener", price = prices[k], variables = mean(ch$n_variables),
               columns = mean(ch$n_columns), cv_auc = a[["estimate"]], cv_auc_se = a[["se"]],
               diff_full = d[["estimate"]], diff_full_se = d[["se"]],
               stability = nogueira_index(freq[, k], M, length(vars)))
  })
  if (include_full) {
    a <- avg(function(r) cv_auc(Y, Pfull[, r], outer[, r]))
    rows[[length(rows) + 1]] <- data.frame(method = "full", price = NA, variables = length(vars), columns = ncol(X),
                                           cv_auc = a[["estimate"]], cv_auc_se = a[["se"]], diff_full = 0,
                                           diff_full_se = 0, stability = NA)
  }
  summ <- do.call(rbind, rows)
  summ$cv_auc_lower <- summ$cv_auc - 1.96 * summ$cv_auc_se
  summ$cv_auc_upper <- summ$cv_auc + 1.96 * summ$cv_auc_se
  summ$diff_full_lower <- summ$diff_full - 1.96 * summ$diff_full_se
  summ$diff_full_upper <- summ$diff_full + 1.96 * summ$diff_full_se
  out <- list(summary = summ, frequency = freq, choices = choices,
              predictions = P, predictions_full = Pfull, folds = outer, prices = prices)
  class(out) <- "nested_superscreener"
  out
}

#' @noRd
.inner_defaults <- function(dots, sl_library, sl_folds) {
  def <- list(K = seq(5, 30, by = 5), rankings = RANKINGS, self_sized = SELF_SIZED, V = 5, loss = "brier",
              library = if (is.null(dots$selection_library)) sl_library else dots$selection_library,
              sl_folds = sl_folds, n_trees = 500)
  dots$selection_library <- NULL
  unknown <- setdiff(names(dots), names(def))
  if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  utils::modifyList(def, dots)
}

#' @export
print.nested_superscreener <- function(x, ...) {
  cat(sprintf("Nested cross-validation of SuperScreener (%d outer folds x %d repeats)\n",
              nrow(unique(x$choices[, c("repeat_id", "fold")])) / ncol(x$folds), ncol(x$folds)))
  s <- x$summary
  show <- data.frame(method = s$method, price = s$price, variables = round(s$variables, 1),
                     cvAUC = sprintf("%.3f (%.3f, %.3f)", s$cv_auc, s$cv_auc_lower, s$cv_auc_upper),
                     diff_from_full = sprintf("%+.3f (%+.3f, %+.3f)", s$diff_full, s$diff_full_lower, s$diff_full_upper),
                     stability = round(s$stability, 2))
  print(show, row.names = FALSE)
  invisible(x)
}
