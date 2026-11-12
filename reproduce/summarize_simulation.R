# summarize_simulation.R
# ---------------------------------------------------------------------------
# Summaries of run_simulation.R: the CSVs behind the chapter's simulation tables and figures,
# and the LaTeX tables themselves.
#
# Fixed budget K (every ranking, and SuperScreener's choice of ranking at that K):
#   pred_all.csv    AUC shortfall from the true risk (x 100), true/noise kept, difference from consensus
#   stab_all.csv    stability index within a replicate (20 training folds) and across replicates
#   choice_all.csv  ranking chosen by SuperScreener at each K (% of replicates)
# Choosing the number of variables (rules applied to the stored risks, no refitting):
#   kselect_picks.csv    one row per replicate x rule x price
#   kselect_summary.csv  means over replicates
#   Rules: ss_full (SuperScreener, every candidate), ss_rank (rankings only), ss_floor (at least
#   5 variables, post hoc), ss_nll (log-likelihood selection), ss_auc (AUC selection), cons_pen
#   (consensus ranking, penalized K), oracle / oracle_auc (selection on the test sample), applied
#   (consensus, smallest K whose 95% lower limit of the CV AUC difference from the full model is
#   above -0.02), sl_screen (screener library, unpenalized), each self-sized screener, full.
# LaTeX: tab_fixedK.tex (Table 2), tab_appendix_rankings.tex (A1), tab_appendix_ksel.tex (A2),
#        tab_appendix_nearempty.tex (A3).
#
# Usage (from the repository root): Rscript reproduce/summarize_simulation.R [--outdir reproduce/output/sim]
# Needs dplyr, tidyr and readr.
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(readr) })
if (requireNamespace("superscreener", quietly = TRUE)) suppressMessages(library(superscreener)) else pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
OUTDIR <- path.expand(arg("--outdir", "reproduce/output/sim"))
SUMDIR <- file.path(OUTDIR, "summary"); dir.create(SUMDIR, showWarnings = FALSE)
C_GRID <- c(0, 0.0025, 0.005, 0.01, 0.02, 0.04)
C_MAIN <- 0.01
K_TAB <- 25L
## Williamson et al. (2025)-style screener library: Spearman and RF rankings at K = 10 and 25,
## screen.glmnet and screen.corP
SL_SCREEN_LIB <- list(spearman = c(10L, 25L), rf = c(10L, 25L), sl_glmnet = 0L, sl_corP = 0L)

reps <- lapply(list.files(file.path(OUTDIR, "reps"), full.names = TRUE), readRDS)
stopifnot(length(reps) > 0)
SCEN <- c(n500 = "Moderate sample", prox_high = "High proxy correlation",
          small_n = "Very small sample", strong_nl = "Strongly nonlinear")
SCEN <- SCEN[names(SCEN) %in% vapply(reps, `[[`, "", "scenario")]
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
pred <- bind_rows(lapply(reps, function(r) mutate(r$pred, scenario = r$scenario, rep = r$rep))) |>
  mutate(short = 100 * (oracle_auc - test_auc))
stab <- bind_rows(lapply(reps, function(r) mutate(r$stab, scenario = r$scenario, rep = r$rep)))
K_RUN <- sort(unique(pred$K[pred$method == "consensus"]))
if (!K_TAB %in% K_RUN) K_TAB <- max(K_RUN)      # e.g. the --quick grid has no K = 25
choice <- bind_rows(lapply(reps, function(r) mutate(r$choice, scenario = r$scenario, rep = r$rep)))

## ---- fixed K -------------------------------------------------------------------------------------
cons <- pred |> filter(method == "consensus") |> select(scenario, rep, K, cons_auc = test_auc)
pred_sum <- pred |> left_join(cons, by = c("scenario", "rep", "K")) |>
  group_by(scenario, method, K) |>
  summarise(n_rep = n(), short_m = mean(short), short_se = mcse(short),
            diff_cons_m = mean(100 * (test_auc - cons_auc)), diff_cons_se = mcse(100 * (test_auc - cons_auc)),
            n_true = mean(n_true), n_noise = mean(n_noise), .groups = "drop")
write_csv(pred_sum, file.path(SUMDIR, "pred_all.csv"))

between <- bind_rows(lapply(split(reps, vapply(reps, `[[`, "", "scenario")), function(rs) {
  p <- rs[[1]]$p
  bind_rows(lapply(names(rs[[1]]$sets), function(k) {
    S <- vapply(rs, function(r) { z <- numeric(p); z[r$sets[[k]]] <- 1; z }, numeric(p))
    tibble(scenario = rs[[1]]$scenario, method = sub(" .*", "", k), K = as.integer(sub(".* ", "", k)),
           phi_between = nogueira_index(rowMeans(matrix(S, nrow = p)), length(rs), p))
  }))
}))
stab_sum <- stab |> group_by(scenario, method, K) |>
  summarise(phi_m = mean(phi, na.rm = TRUE), phi_q1 = quantile(phi, 0.25, na.rm = TRUE),
            phi_q3 = quantile(phi, 0.75, na.rm = TRUE), stable_n = mean(stable_n), stable_true = mean(stable_true),
            stable_noise = mean(stable_noise), .groups = "drop") |>
  left_join(between, by = c("scenario", "method", "K"))
write_csv(stab_sum, file.path(SUMDIR, "stab_all.csv"))
choice_sum <- choice |> count(scenario, K, chosen) |> group_by(scenario, K) |>
  mutate(pct = 100 * n / sum(n)) |> ungroup()
write_csv(choice_sum, file.path(SUMDIR, "choice_all.csv"))

## ---- choosing the number of variables ----------------------------------------------------------
pick_rep <- function(r) {
  pr <- r$pred
  cand <- pr |> filter((method %in% RANKINGS & K > 0) | method %in% SELF_SIZED) |>
    transmute(method, K = ifelse(method %in% SELF_SIZED, NA_integer_, K), size = n_cols, cv_risk = cv_brier,
              cv_nll, cv_auc, test_risk = test_brier, test_auc, oracle_auc, n_true, n_noise)
  ## the oracle applies the rule to the test-sample risks, with the price set from the test risk
  ## of the intercept-only model
  oracle_at <- function(c) choose_candidate(mutate(cand, cv_risk = test_risk), c, r$test_null_brier)
  pick_auc <- function(cd, c, col) { g <- max(cd[[col]]) - 0.5; cd[order(-(cd[[col]] - c * g * cd$size), cd$size)[1], ] }
  row_out <- function(rule, c, row, ref) {
    tibble(rule = rule, c = c, method = row$method, K = row$K, n_vars = row$size,
           short = 100 * (row$oracle_auc - row$test_auc), n_true = row$n_true, n_noise = row$n_noise,
           regret = 100 * (row$oracle_auc - row$test_auc) - 100 * (ref$oracle_auc - ref$test_auc))
  }
  out <- list()
  for (c in C_GRID) {
    o <- oracle_at(c); oa <- pick_auc(cand, c, "test_auc")
    out <- c(out, list(
      row_out("ss_full", c, choose_candidate(cand, c, r$cv_null_brier), o),
      row_out("ss_rank", c, choose_candidate(filter(cand, method %in% RANKINGS), c, r$cv_null_brier), o),
      row_out("ss_floor", c, choose_candidate(cand, c, r$cv_null_brier, min_size = 5), o),
      row_out("cons_pen", c, choose_candidate(filter(cand, method == "consensus"), c, r$cv_null_brier), o),
      row_out("ss_nll", c, choose_candidate(mutate(cand, cv_risk = cv_nll), c, r$cv_null_nll), o),
      row_out("oracle", c, o, o),
      row_out("ss_auc", c, pick_auc(cand, c, "cv_auc"), oa),
      row_out("oracle_auc", c, oa, o)))
  }
  consK <- pr |> filter(method == "consensus", K > 0) |> arrange(K)
  ok <- consK$cv_diff_full - 1.96 * consK$cv_diff_full_se > -0.02
  kA <- if (any(ok)) consK$K[which(ok)[1]] else max(consK$K)
  in_lib <- mapply(function(m, k) !is.null(SL_SCREEN_LIB[[m]]) && (if (is.na(k)) 0L else k) %in% SL_SCREEN_LIB[[m]],
                   cand$method, cand$K)
  fixed <- list(applied = filter(cand, method == "consensus", K == kA),
                sl_screen = cand[in_lib, ] |> slice_min(cv_risk, n = 1, with_ties = FALSE))
  for (m in SELF_SIZED) fixed[[m]] <- filter(cand, method == m)
  fu <- filter(pr, method == "full")
  fixed$full <- tibble(method = "full", K = NA_integer_, size = fu$n_cols, test_auc = fu$test_auc,
                       oracle_auc = fu$oracle_auc, n_true = fu$n_true, n_noise = fu$n_noise)
  for (nm in names(fixed)) for (c in C_GRID) out[[length(out) + 1]] <- row_out(nm, c, fixed[[nm]], oracle_at(c))
  st <- r$stab |> transmute(method, K = ifelse(method %in% SELF_SIZED, NA_integer_, K), phi)
  bind_rows(out) |> left_join(st, by = c("method", "K")) |> mutate(scenario = r$scenario, rep = r$rep)
}
picks <- bind_rows(lapply(reps, pick_rep))
write_csv(picks, file.path(SUMDIR, "kselect_picks.csv"))
ksum <- picks |> group_by(scenario, rule, c) |>
  summarise(n_rep = n(), vars_m = mean(n_vars), vars_sd = sd(n_vars), short_m = mean(short), short_se = mcse(short),
            true_m = mean(n_true), noise_m = mean(n_noise), regret_m = mean(regret), regret_se = mcse(regret),
            phi_m = mean(phi, na.rm = TRUE), .groups = "drop")
write_csv(ksum, file.path(SUMDIR, "kselect_summary.csv"))

## ---- LaTeX tables ------------------------------------------------------------------------------------
fmt <- function(x, f) ifelse(is.na(x) | is.nan(x), "", sprintf(f, x))
SCEN_TEX <- c(n500 = "Moderate sample ($n=500$, $p=300$)", prox_high = "High proxy correlation ($n=1000$, $p=300$)",
              small_n = "Very small sample ($n=100$, $p=2000$)", strong_nl = "Strongly nonlinear ($n=200$, $p=500$)")
RK <- c(lrt = "Univariate LRT", spearman = "Spearman (\\texttt{screen.corRank})", rf = "Random forest importance",
        lasso = "Lasso entry order", consensus = "Consensus ranking")

## Table 2: AUC shortfall and stability [IQR] at K = 25
ROWS2 <- c(consensus = "Consensus ranking (applied analysis)", ss = "SuperScreener (best ranking at $K$)", RK[c("lrt", "spearman", "rf", "lasso")],
           full = "Full model (all columns)")
cell2 <- function(sc, m) {
  sh <- pred_sum$short_m[pred_sum$scenario == sc & pred_sum$method == m & (pred_sum$K %in% K_TAB | m == "full")]
  st <- stab_sum[stab_sum$scenario == sc & stab_sum$method == m & stab_sum$K == K_TAB, ]
  c(fmt(sh, "%.2f"), if (nrow(st)) sprintf("%.2f [%.2f, %.2f]", st$phi_m, st$phi_q1, st$phi_q3) else "")
}
writeLines(c(sprintf("\\begin{tabular}{@{}l*{%d}{rr}@{}}", length(SCEN)), "\\toprule",
             paste0(" & ", paste(sprintf("\\multicolumn{2}{c}{%s}", SCEN), collapse = " & "), " \\\\"),
             paste0("Method", strrep(" & Shortfall & Stability", length(SCEN)), " \\\\"), "\\midrule",
             vapply(names(ROWS2), function(m) paste0(paste(c(ROWS2[[m]], unlist(lapply(names(SCEN), cell2, m = m))), collapse = " & "), " \\\\"), ""),
             "\\bottomrule", "\\end{tabular}"), file.path(SUMDIR, "tab_fixedK.tex"))

## Table A1: each ranking at K = 25
rowsA1 <- unlist(lapply(names(SCEN), function(sc) {
  c(sprintf("\\multicolumn{6}{@{}l}{\\textit{%s}} \\\\", SCEN_TEX[[sc]]), vapply(names(RK), function(m) {
    pr <- pred_sum[pred_sum$scenario == sc & pred_sum$method == m & pred_sum$K == K_TAB, ]
    st <- stab_sum[stab_sum$scenario == sc & stab_sum$method == m & stab_sum$K == K_TAB, ]
    pc <- choice_sum$pct[choice_sum$scenario == sc & choice_sum$K == K_TAB & choice_sum$chosen == m]
    sprintf("\\quad %s & %.0f & %.1f & %.1f & %.2f & %.2f \\\\", RK[[m]], if (length(pc)) pc else 0,
            pr$n_true, pr$n_noise, st$phi_m, st$phi_between)
  }, ""))
}))
writeLines(c("\\begin{tabular}{@{}lrrrrr@{}}", "\\toprule",
             " & Chosen by & True & Noise & \\multicolumn{2}{c}{Stability} \\\\", "\\cmidrule(l){5-6}",
             "Ranking & SuperScreener (\\%) & kept & kept & within & across \\\\", "\\midrule",
             rowsA1, "\\bottomrule", "\\end{tabular}"), file.path(SUMDIR, "tab_appendix_rankings.tex"))

## Table A2: choosing the number of variables at c = 1%
RULES <- list(c("ss_full", C_MAIN, "SuperScreener, pre-specified"), c("ss_floor", C_MAIN, "SuperScreener, at least 5 variables (post hoc)"),
              c("ss_rank", C_MAIN, "SuperScreener, rankings only"), c("ss_auc", C_MAIN, "SuperScreener, AUC selection"),
              c("ss_full", 0, "SuperScreener, no price ($c=0$)"), c("applied", C_MAIN, "Rule of the applied analysis"),
              c("sl_screen", C_MAIN, "Screener library (unpenalized)"), c("sl_corP", C_MAIN, "\\texttt{screen.corP}"),
              c("sl_glmnet", C_MAIN, "\\texttt{screen.glmnet}"), c("full", C_MAIN, "Full model (all columns)"),
              c("oracle", C_MAIN, "Oracle (same price)"))
rowsA2 <- unlist(lapply(names(SCEN), function(sc) {
  c(sprintf("\\multicolumn{6}{@{}l}{\\textit{%s}} \\\\", SCEN_TEX[[sc]]), vapply(RULES, function(r) {
    x <- ksum[ksum$scenario == sc & ksum$rule == r[1] & abs(ksum$c - as.numeric(r[2])) < 1e-9, ]
    o <- ksum[ksum$scenario == sc & ksum$rule == "oracle" & abs(ksum$c - as.numeric(r[2])) < 1e-9, ]
    sprintf("\\quad %s & %.1f & %.1f & %.2f & %.2f & %s \\\\", r[3], x$vars_m, x$vars_sd, x$short_m,
            x$short_m - o$short_m, if (r[1] == "full") "" else fmt(x$phi_m, "%.2f"))
  }, ""))
}))
writeLines(c("\\begin{tabular}{@{}lrrrrr@{}}", "\\toprule",
             " & \\multicolumn{2}{c}{Variables used} & AUC & Gap to & \\\\", "\\cmidrule(lr){2-3}",
             "Rule & Mean & SD & shortfall & oracle & Stability \\\\", "\\midrule",
             rowsA2, "\\bottomrule", "\\end{tabular}"), file.path(SUMDIR, "tab_appendix_ksel.tex"))

## Table A3: near-empty candidates
per_rep <- bind_rows(lapply(reps, function(r) {
  pr <- r$pred
  cand <- pr |> filter((method %in% RANKINGS & K > 0) | method %in% SELF_SIZED)
  ss <- function(m) pr$n_cols[pr$method == m]
  tibble(scenario = r$scenario, none_beats = all(cand$cv_brier[cand$n_cols >= 1] >= r$cv_null_brier),
         l1_zero = ss("lasso_1se") == 0, l1_n = ss("lasso_1se"), sb_zero = ss("stabsel") == 0,
         sb_n = ss("stabsel"), gl_lt5 = ss("sl_glmnet") < 5)
}))
pk <- picks |> filter(rule == "ss_full", abs(c - C_MAIN) < 1e-9) |> group_by(scenario) |>
  summarise(ch0 = 100 * mean(n_vars == 0), ch5 = 100 * mean(n_vars < 5),
            sh5 = if (any(n_vars < 5)) mean(short[n_vars < 5]) else NA_real_)
sm <- per_rep |> group_by(scenario) |>
  summarise(none = 100 * mean(none_beats), l1z = 100 * mean(l1_zero), l1n = mean(l1_n), sbz = 100 * mean(sb_zero),
            sbn = mean(sb_n), gl = 100 * mean(gl_lt5)) |> left_join(pk, by = "scenario")
sm <- sm[match(names(SCEN), sm$scenario), ]
line3 <- function(lab, v, f) paste0(lab, " & ", paste(fmt(v, f), collapse = " & "), " \\\\")
writeLines(c(sprintf("\\begin{tabular}{@{}l%s@{}}", strrep("r", length(SCEN))), "\\toprule",
             paste0(" & ", paste(SCEN, collapse = " & "), " \\\\"), "\\midrule",
             line3("No candidate with a variable beats the intercept-only model (\\%)", sm$none, "%.0f"), "\\addlinespace",
             line3("Lasso at \\texttt{lambda.1se} keeps no variables (\\%)", sm$l1z, "%.0f"),
             line3("Lasso at \\texttt{lambda.1se}, mean variables kept", sm$l1n, "%.1f"),
             line3("Stability selection keeps no variables (\\%)", sm$sbz, "%.0f"),
             line3("Stability selection, mean variables kept", sm$sbn, "%.1f"),
             line3("\\texttt{screen.glmnet} keeps two to four variables (\\%)", sm$gl, "%.0f"), "\\addlinespace",
             line3("SuperScreener at $c=1\\%$ chooses no variables (\\%)", sm$ch0, "%.0f"),
             line3("SuperScreener at $c=1\\%$ chooses fewer than five (\\%)", sm$ch5, "%.0f"),
             line3("\\quad AUC shortfall in those datasets", sm$sh5, "%.1f"),
             "\\bottomrule", "\\end{tabular}"), file.path(SUMDIR, "tab_appendix_nearempty.tex"))

message(sprintf("%d replicates summarized -> %s", length(reps), SUMDIR))
print(table(vapply(reps, `[[`, "", "scenario")))
