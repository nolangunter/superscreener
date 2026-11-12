# apply_to_your_data.R
# ---------------------------------------------------------------------------
# Template for applying SuperScreener to your own data, as in the chapter's registry case study:
# the price path (how many variables are kept at each price c), nested cross-validated AUC of
# the whole procedure with 95% intervals, the difference from the model with every column, and
# how often each variable is selected across outer training folds.
#
# Your data never leave your machine, and every file written below holds aggregates only.
#
#   Rscript examples/apply_to_your_data.R --data path/to/data.csv --outcome died_30d \
#       [--drop id,date] [--outdir my_results] [--cores 4] [--repeats 1] [--fast]
#
# --fast is a trial run (light learners: glm, lasso, mean) to check that everything works
# before the real run.
#
# Requirements for the CSV: one row per subject, a 0/1 outcome column, predictors that are
# numeric, logical, character or factor. Missing values must be imputed first. Categorical
# variables are expanded to indicator columns and charged as ONE variable each.
#
# Run time: each outer fold reruns every candidate (34 by default) in each of 5 inner folds,
# so the nested step costs roughly 5 x 5 x 34 Super Learner fits per repeat. On a laptop,
# a few hundred rows and a few dozen variables take on the order of an hour with --cores 5;
# registry-sized data (thousands of rows, hundreds of columns) belong on a server.
# ---------------------------------------------------------------------------

## one thread per process, so parallel outer folds do not compete for cores
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
library(superscreener)
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
DATA <- arg("--data", NA); OUTCOME <- arg("--outcome", NA)
if (is.na(DATA) || is.na(OUTCOME)) stop("give --data <csv> and --outcome <column>")
DROP <- strsplit(arg("--drop", ""), ",")[[1]]
OUTDIR <- arg("--outdir", "superscreener_results"); dir.create(OUTDIR, showWarnings = FALSE)
CORES <- as.integer(arg("--cores", "1")); REPEATS <- as.integer(arg("--repeats", "1"))
PRICES <- c(0, 0.0025, 0.005, 0.01, 0.02, 0.04)

dat <- utils::read.csv(DATA, stringsAsFactors = TRUE, check.names = FALSE)
Y <- dat[[OUTCOME]]
pred_df <- dat[, setdiff(names(dat), c(OUTCOME, DROP)), drop = FALSE]
ex <- expand_factors(pred_df)                  # numeric columns + column -> variable map
message(sprintf("%d rows, %d variables, %d columns, %.0f%% cases", nrow(ex$X), length(unique(ex$groups)),
                ncol(ex$X), 100 * mean(Y)))

## Inner selection library: a smaller learner keeps the nested computation feasible
## (the registry analysis used glmnet, ranger and the mean). The final models use sl_library.
FAST <- "--fast" %in% args
SEL_LIB <- if (FAST) c("SL.glm", "SL.glmnet", "SL.mean") else c("SL.glmnet", "SL.ranger_300", "SL.mean")
FINAL_LIB <- if (FAST) SEL_LIB else sl_library_default()

## 1. SuperScreener on all the data: the recommended set at each price
fit <- superscreener(ex$X, Y, price = 0.01, groups = ex$groups,
                     selection_library = SEL_LIB, sl_library = FINAL_LIB, seed = 1)
print(fit)
pp <- price_path(fit, PRICES)
pp$variables <- vapply(pp$candidate, function(k) paste(unique(ex$groups[fit$sets[[k]]]), collapse = "; "), "")
utils::write.csv(pp, file.path(OUTDIR, "price_path.csv"), row.names = FALSE)

## 2. Nested cross-validation of the whole procedure (the honest accuracy)
ncv <- nested_cv_superscreener(ex$X, Y, prices = PRICES, groups = ex$groups, V_outer = 5, repeats = REPEATS,
                               selection_library = SEL_LIB, sl_library = FINAL_LIB, cores = CORES, seed = 1)
print(ncv)
utils::write.csv(ncv$summary, file.path(OUTDIR, "nested_cv_summary.csv"), row.names = FALSE)
freq <- data.frame(variable = rownames(ncv$frequency), ncv$frequency, check.names = FALSE)
utils::write.csv(freq[order(-freq[["c=0.01"]]), ], file.path(OUTDIR, "selection_frequency.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(table(price = ncv$choices$price, candidate = ncv$choices$candidate)),
                 file.path(OUTDIR, "candidates_chosen.csv"), row.names = FALSE)

## 3. Price path figure (needs ggplot2)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  s <- ncv$summary[ncv$summary$method == "superscreener" & ncv$summary$price > 0, ]
  g <- ggplot(s, aes(100 * price, 100 * diff_full)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = 100 * diff_full_lower, ymax = 100 * diff_full_upper), width = 0, color = "#eb6834") +
    geom_line(color = "#eb6834") + geom_point(color = "#eb6834") +
    geom_text(aes(label = sprintf("%.1f vars", variables)), vjust = -1, size = 3) +
    scale_x_log10() + theme_minimal() +
    labs(x = "Price c (% of the best candidate's risk reduction, log scale)",
         y = "Nested cvAUC minus full model (x 100)")
  ggsave(file.path(OUTDIR, "price_path.png"), g, width = 6, height = 4, dpi = 200, bg = "white")
}
message("results (aggregates only) in ", OUTDIR)
