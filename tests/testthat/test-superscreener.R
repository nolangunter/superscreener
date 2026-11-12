d <- make_design(p = 40, mc_n = 2e4)
dat <- simulate_data(250, d, seed = 2)
lib <- c("SL.glm", "SL.mean")

test_that("superscreener runs end to end and is reproducible", {
  f1 <- superscreener(dat$X, dat$Y, K = c(5, 10), sl_library = lib, seed = 1, verbose = FALSE)
  f2 <- superscreener(dat$X, dat$Y, K = c(5, 10), sl_library = lib, seed = 1, verbose = FALSE)
  expect_s3_class(f1, "superscreener")
  expect_equal(nrow(f1$candidates), length(RANKINGS) * 2 + length(SELF_SIZED))
  expect_equal(f1$candidates$cv_risk, f2$candidates$cv_risk)
  expect_equal(f1$selected, f1$sets[[f1$chosen$candidate]])
  # the SuperLearner fits did not silently fall back to the intercept
  expect_gt(sd(f1$cv_predictions[, "lrt_K5"]), 0.01)
  p <- predict(f1, dat$X)
  expect_length(p, 250); expect_true(all(p > 0 & p < 1))
  pp <- price_path(f1, c(0, 0.01, 0.1))
  expect_true(all(diff(pp$size) <= 0))
  f3 <- set_price(f1, 0.1)
  expect_equal(f3$chosen$candidate, pp$candidate[3])
})

test_that("groups make the cost count variables, not columns", {
  g <- setNames(rep(paste0("v", 1:20), each = 2), colnames(dat$X))
  f <- superscreener(dat$X, dat$Y, K = 10, rankings = "lrt", self_sized = character(0),
                     sl_library = lib, groups = g, seed = 1, verbose = FALSE)
  expect_equal(f$candidates$size, length(unique(g[f$sets$lrt_K10])))
})

test_that("input checks", {
  expect_error(superscreener(dat$X, dat$Y + 1, verbose = FALSE), "0/1")
  X <- dat$X; X$chr <- "a"
  expect_error(superscreener(X, dat$Y, verbose = FALSE), "expand_factors")
})

test_that("nested cross-validation returns a summary per price", {
  ncv <- nested_cv_superscreener(dat$X, dat$Y, prices = c(0, 0.05), K = c(5, 10), rankings = c("lrt", "spearman"),
                                 self_sized = "sl_corP", sl_library = lib, V_outer = 3, seed = 1, verbose = FALSE)
  expect_s3_class(ncv, "nested_superscreener")
  expect_equal(nrow(ncv$summary), 3)
  expect_true(all(ncv$summary$cv_auc > 0.5))
  expect_equal(dim(ncv$frequency), c(40, 2))
})

test_that("a fit where only SL.mean survives predicts the case proportion, not 0", {
  p <- fit_sl(dat$X[, "true1", drop = FALSE], dat$Y, dat$X[1:3, "true1", drop = FALSE], c("SL.glmnet", "SL.mean"))
  expect_equal(p, rep(mean(dat$Y), 3))
})

test_that("parallel nested CV gives the same answer and handles arguments passed as variables", {
  skip_on_cran()
  run <- function(cores) {
    my_prices <- c(0, 0.05)
    nested_cv_superscreener(dat$X, dat$Y, prices = my_prices, K = 5, rankings = "lrt", self_sized = "sl_corP",
                            sl_library = lib, V_outer = 2, seed = 3, cores = cores, verbose = FALSE)$summary
  }
  expect_equal(run(2), run(1))
})
