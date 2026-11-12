# run_simulation.R
# ---------------------------------------------------------------------------
# Simulation study of the SuperScreener chapter: one file per (scenario, replicate).
#
# Per replicate (training sample of size n, independent test sample of 10,000):
#   1. Every ranking on the whole training sample and on each of 4 x 5 cross-validation
#      training folds (4 repeated 5-fold splits); every self-sized screener likewise.
#   2. Cross-validation on the first split: the Super Learner fit to each candidate's set
#      in each training fold, predicting the held-out fold. Gives each candidate's
#      cross-validated Brier risk, log-likelihood and AUC (and the AUC difference from the
#      model with every column, used by the applied rule).
#   3. Every candidate, the full model and the true set fit to the whole training sample
#      and scored on the test sample. "ss" at budget K is the ranking with the smallest
#      cross-validated Brier risk at that K.
#   4. Stability of each candidate's set across the 20 training folds (Nogueira index).
# Choosing K (SuperScreener's penalized rule, the applied rule, the screener library and
# the oracle) needs no further fitting: summarize_simulation.R applies the rules to the
# stored cross-validated and test risks.
#
# Usage (from the repository root, after installing the package):
#   Rscript reproduce/run_simulation.R --quick                    # smoke test: small designs, light learner
#   Rscript reproduce/run_simulation.R --reps 1000 --cores 50     # the chapter run
#     [--scenarios n500,prox_high,small_n,strong_nl] [--outdir reproduce/output/sim]
# Resumable: a replicate whose file exists is skipped, and a larger --reps adds replicates.
# Cost per replicate on one core (2026 server): n500 ~10 min, prox_high ~40, small_n ~5,
# strong_nl ~4, so 1000 replicates of all four is about 1000 core-hours.
# ---------------------------------------------------------------------------

if (requireNamespace("superscreener", quietly = TRUE)) {
  suppressMessages(library(superscreener))
} else {
  pkgload::load_all(".", quiet = TRUE)     # run from the repository root without installing
}
suppressMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
QUICK <- "--quick" %in% args
CORES <- as.integer(arg("--cores", if (QUICK) "4" else "1"))
REPS <- seq_len(as.integer(arg("--reps", if (QUICK) "2" else "10")))
OUTDIR <- path.expand(arg("--outdir", file.path("reproduce", "output", if (QUICK) "sim_quick" else "sim")))
R_STAB <- if (QUICK) 2L else 4L          # repeated 5-fold splits for stability
V <- 5L
K_GRID <- if (QUICK) c(5L, 10L, 20L) else c(5L, 10L, 15L, 20L, 25L, 30L)
N_TEST <- if (QUICK) 2000L else 10000L
V_INNER <- if (QUICK) 3L else 5L         # the Super Learner's own folds
SL_LIB <- if (QUICK) c("SL.glm", "SL.glmnet", "SL.mean") else sl_library_default()
PI_THR <- 0.8                             # stable set: selected in >= 80% of training folds
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
dir.create(file.path(OUTDIR, "reps"), recursive = TRUE, showWarnings = FALSE)

SCENARIOS <- chapter_scenarios()
if (QUICK) SCENARIOS <- list(
  prox_high = list(n = 400, p = 60, prox_cor = c(0.90, 0.97), surface = "nonlinear"),
  strong_nl = list(n = 200, p = 80, prox_cor = c(0.5, 0.85), surface = "strong_nonlinear"))
SCEN <- strsplit(arg("--scenarios", paste(names(SCENARIOS), collapse = ",")), ",")[[1]]
stopifnot(all(SCEN %in% names(SCENARIOS)))
## seeds depend on the scenario's position in the chapter's original list, which had two more
## scenarios ("base", "linear") after these four
SCEN_INDEX <- c(n500 = 1, prox_high = 2, small_n = 3, strong_nl = 4)

designs_file <- file.path(OUTDIR, "designs.rds")
designs <- if (file.exists(designs_file)) readRDS(designs_file) else list()
for (nm in setdiff(SCEN, names(designs))) {
  message("calibrating design ", nm)
  s <- SCENARIOS[[nm]]
  designs[[nm]] <- make_design(p = s$p, prox_cor = s$prox_cor, surface = s$surface, name = nm)
  designs[[nm]]$cell <- s
}
saveRDS(designs, designs_file)

roles_count <- function(cols, roles) {
  r <- roles[cols]
  c(true = sum(r == "true"), proxy = sum(r == "proxy"), noise = sum(r == "noise"))
}

stab_row <- function(fr, method, K, p, d, full_set, M) {
  stable <- names(fr)[fr >= PI_THR]
  unstable <- intersect(full_set, names(fr)[fr < 0.5])
  rs <- roles_count(stable, d$roles); ru <- roles_count(unstable, d$roles)
  data.frame(method = method, K = K, phi = nogueira_index(fr, M, p), n_distinct = sum(fr > 0),
             stable_n = length(stable), stable_true = rs[["true"]], stable_proxy = rs[["proxy"]],
             stable_noise = rs[["noise"]], unstable_n = length(unstable),
             unstable_proxy = ru[["proxy"]], unstable_noise = ru[["noise"]], freq_true = mean(fr[d$S0]))
}

run_rep <- function(scen, rep) {
  f_out <- file.path(OUTDIR, "reps", sprintf("%s_r%03d.rds", scen, rep))
  if (file.exists(f_out)) return("skip")
  t0 <- Sys.time()
  d <- designs[[scen]]
  seed <- 1e5 * SCEN_INDEX[[scen]] + rep
  tr <- simulate_data(d$cell$n, d, seed)
  te <- simulate_data(N_TEST, d, seed + 5e4)
  X <- tr$X; Y <- tr$Y; cols <- colnames(X); p <- length(cols)
  set.seed(seed + 1)
  folds <- sapply(seq_len(R_STAB), function(r) stratified_folds(Y, V))
  fold_id <- function(r, v) sprintf("r%d_f%d", r, v)

  ## 1. rankings and self-sized sets: whole training sample and every training fold
  rk_full <- rank_columns(X, Y)
  rk_fold <- list()
  for (r in seq_len(R_STAB)) for (v in seq_len(V)) {
    trn <- folds[, r] != v
    rk_fold[[fold_id(r, v)]] <- rank_columns(X[trn, , drop = FALSE], Y[trn])
  }
  sz_full <- self_sized_sets(X, Y, seed = seed + 2)
  sz_fold <- list()
  for (r in seq_len(R_STAB)) for (v in seq_len(V)) {
    trn <- folds[, r] != v
    sz_fold[[fold_id(r, v)]] <- self_sized_sets(X[trn, , drop = FALSE], Y[trn], seed = seed + 10 * r + v)
  }

  ## 2. cross-validation on the first split (K = 0 marks a self-sized screener)
  f1 <- folds[, 1]
  oof <- list(); oof_full <- oof_null <- numeric(length(Y))
  for (v in seq_len(V)) {
    trn <- f1 != v; tst <- !trn
    rk <- rk_fold[[fold_id(1, v)]]
    oof_full[tst] <- fit_sl(X[trn, , drop = FALSE], Y[trn], X[tst, , drop = FALSE], SL_LIB, V_INNER)
    oof_null[tst] <- mean(Y[trn])
    for (m in RANKINGS) for (K in K_GRID) {
      key <- paste(m, K); if (is.null(oof[[key]])) oof[[key]] <- numeric(length(Y))
      cc <- top_k(rk[, m], K)
      oof[[key]][tst] <- fit_sl(X[trn, cc, drop = FALSE], Y[trn], X[tst, cc, drop = FALSE], SL_LIB, V_INNER)
    }
    for (m in SELF_SIZED) {
      key <- paste(m, 0L); if (is.null(oof[[key]])) oof[[key]] <- numeric(length(Y))
      cc <- sz_fold[[fold_id(1, v)]][[m]]
      oof[[key]][tst] <- fit_sl(X[trn, cc, drop = FALSE], Y[trn], X[tst, cc, drop = FALSE], SL_LIB, V_INNER)
    }
  }
  cv_row <- function(method, K, pr) {
    a <- cv_auc(Y, pr, f1); dfull <- cv_auc(Y, pr, f1, p_ref = oof_full)
    data.frame(method = method, K = K, cv_brier = brier(Y, pr), cv_nll = neg_loglik(Y, pr),
               cv_auc = a$estimate, cv_auc_se = a$se, cv_diff_full = dfull$estimate, cv_diff_full_se = dfull$se)
  }
  cv_tab <- do.call(rbind, lapply(names(oof), function(key)
    cv_row(sub(" .*", "", key), as.integer(sub(".* ", "", key)), oof[[key]])))

  ## 3. whole-sample fits scored on the test sample
  test_row <- function(method, K, cc) {
    pr <- fit_sl(X[, cc, drop = FALSE], Y, te$X[, cc, drop = FALSE], SL_LIB, V_INNER)
    rc <- roles_count(cc, d$roles)
    data.frame(method = method, K = K, n_cols = length(cc), n_true = rc[["true"]], n_proxy = rc[["proxy"]],
               n_noise = rc[["noise"]], test_auc = auc(te$Y, pr), test_nll = neg_loglik(te$Y, pr),
               test_brier = brier(te$Y, pr))
  }
  sets <- list(); pred <- list(); choice <- list()
  for (K in K_GRID) {
    cvK <- cv_tab[cv_tab$K == K & cv_tab$method %in% RANKINGS, ]
    chosen <- cvK$method[which.min(cvK$cv_brier)]
    choice[[length(choice) + 1]] <- data.frame(K = K, chosen = chosen)
    for (m in RANKINGS) {
      cc <- top_k(rk_full[, m], K)
      pred[[length(pred) + 1]] <- test_row(m, K, cc); sets[[paste(m, K)]] <- cc
    }
    ss <- pred[[which(vapply(pred, function(z) z$method == chosen && z$K == K, logical(1)))]]
    ss$method <- "ss"; pred[[length(pred) + 1]] <- ss; sets[[paste("ss", K)]] <- sets[[paste(chosen, K)]]
  }
  for (m in SELF_SIZED) { pred[[length(pred) + 1]] <- test_row(m, 0L, sz_full[[m]]); sets[[paste(m, 0L)]] <- sz_full[[m]] }
  pred <- do.call(rbind, c(pred, list(test_row("full", NA_integer_, cols), test_row("true_set", length(d$S0), d$S0))))
  choice <- do.call(rbind, choice)
  cv_ss <- merge(choice, cv_tab, by.x = c("chosen", "K"), by.y = c("method", "K"))
  cv_ss$method <- "ss"; cv_ss$chosen <- NULL
  cv_all <- rbind(cv_tab, cv_ss[, names(cv_tab)], cv_row("full", NA_integer_, oof_full))
  pred <- merge(pred, cv_all, by = c("method", "K"), all.x = TRUE, sort = FALSE)
  pred$oracle_auc <- auc(te$Y, te$p)

  ## 4. stability across the R_STAB x V training folds
  M <- length(rk_fold)
  stab <- list()
  ## As in the chapter, a fold's set is every column with rank <= K. This differs from the top K
  ## only for the lasso ranking, whose entry ranks can tie.
  for (m in RANKINGS) for (K in K_GRID) {
    fr <- selection_frequency(lapply(rk_fold, function(z) names(which(z[, m] <= K))), cols)
    stab[[length(stab) + 1]] <- stab_row(fr, m, K, p, d, top_k(rk_full[, m], K), M)
  }
  for (m in SELF_SIZED) {
    fr <- selection_frequency(lapply(sz_fold, `[[`, m), cols)
    stab[[length(stab) + 1]] <- stab_row(fr, m, 0L, p, d, sz_full[[m]], M)
  }
  stab <- do.call(rbind, stab)
  ss_stab <- merge(choice, stab, by.x = c("chosen", "K"), by.y = c("method", "K"))
  ss_stab$method <- "ss"; ss_stab$chosen <- NULL
  stab <- rbind(stab, ss_stab[, names(stab)])

  res <- list(scenario = scen, rep = rep, seed = seed, pred = pred, stab = stab, choice = choice,
              cv_null_brier = brier(Y, oof_null), cv_null_nll = neg_loglik(Y, oof_null),
              test_null_brier = brier(te$Y, rep(mean(Y), length(te$Y))),
              sets = lapply(sets, function(s) match(s, cols)), n = d$cell$n, p = p,
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(res, f_out)
  sprintf("%s rep %d: %.1f min", scen, rep, res$minutes)
}

tasks <- expand.grid(rep = REPS, scen = SCEN, stringsAsFactors = FALSE)
tasks <- tasks[order(tasks$rep, match(tasks$scen, SCEN)), ]
message(sprintf("%d replicates on %d cores -> %s", nrow(tasks), CORES, OUTDIR))
out <- mclapply(seq_len(nrow(tasks)), function(i) {
  tryCatch(run_rep(tasks$scen[i], tasks$rep[i]), error = function(e) {
    msg <- sprintf("%s rep %d FAILED: %s", tasks$scen[i], tasks$rep[i], conditionMessage(e))
    cat(msg, "\n", file = file.path(OUTDIR, "errors.log"), append = TRUE); msg
  })
}, mc.cores = CORES, mc.preschedule = FALSE)
message(paste(unlist(out), collapse = "\n"))
