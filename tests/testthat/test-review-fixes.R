## Regression tests for issues found in the pre-release review.
d <- make_design(p = 40, mc_n = 2e4)
dat <- simulate_data(250, d, seed = 2)
lib <- c("SL.glm", "SL.mean")

test_that("nested CV survives an empty chosen set", {
  set.seed(9)
  X <- data.frame(matrix(rnorm(150 * 20), 150)); Y <- rbinom(150, 1, 0.3)   # pure noise
  ncv <- nested_cv_superscreener(X, Y, prices = c(0, 5), K = 5, rankings = "lrt", self_sized = "lasso_1se",
                                 sl_library = lib, V_outer = 3, seed = 1, verbose = FALSE)
  expect_false(anyNA(ncv$predictions))
})

test_that("factor groups charge the right cost", {
  X <- dat$X[, 1:4]
  g <- factor(c("b", "b", "a", "c")); names(g) <- colnames(X)
  f <- superscreener(X, dat$Y, K = 1, rankings = "lrt", self_sized = character(0), sl_library = lib,
                     groups = g, cost = c(a = 1, b = 100, c = 1), seed = 1, verbose = FALSE)
  first <- f$sets$lrt_K1
  expect_equal(f$candidates$size, unname(c(a = 1, b = 100, c = 1)[as.character(g[first])]))
})

test_that("logical and factor outcomes are converted", {
  f1 <- superscreener(dat$X, dat$Y, K = 5, rankings = "lrt", self_sized = character(0), sl_library = lib,
                      seed = 1, verbose = FALSE)
  f2 <- superscreener(dat$X, dat$Y == 1, K = 5, rankings = "lrt", self_sized = character(0), sl_library = lib,
                      seed = 1, verbose = FALSE)
  f3 <- superscreener(dat$X, factor(dat$Y), K = 5, rankings = "lrt", self_sized = character(0), sl_library = lib,
                      seed = 1, verbose = FALSE)
  expect_equal(f2$candidates$cv_risk, f1$candidates$cv_risk)
  expect_equal(f3$candidates$cv_risk, f1$candidates$cv_risk)
  expect_error(superscreener(dat$X, factor(ifelse(dat$Y == 1, "yes", "no")), verbose = FALSE), "0/1")
})

test_that("expand_factors refuses to overwrite a column", {
  df <- data.frame(sex_m = c(1, 2, 3, 4), sex = c("f", "m", "f", "m"))
  expect_error(expand_factors(df), "collide")
})

test_that("seed = NULL works in nested CV", {
  expect_silent(nested_cv_superscreener(dat$X, dat$Y, prices = 0, K = 5, rankings = "lrt", self_sized = character(0),
                                        sl_library = lib, V_outer = 2, seed = NULL, verbose = FALSE))
})

test_that("too few cases for the folds is an error, not NaN", {
  Y <- c(rep(1, 3), rep(0, 47)); X <- dat$X[1:50, 1:5]
  expect_error(superscreener(X, Y, verbose = FALSE), "fewer than the 5 folds")
})

test_that("budgets above p are capped and deduplicated", {
  X <- dat$X[, 1:8]
  f <- superscreener(X, dat$Y, K = c(5, 10, 20), rankings = "lrt", self_sized = character(0), sl_library = lib,
                     seed = 1, verbose = FALSE)
  expect_equal(f$candidates$K, c(5L, 8L))
})

test_that("bad arguments give clear errors", {
  expect_error(superscreener(dat$X, dat$Y, rankings = "las", verbose = FALSE), "Unknown ranking")
  expect_error(superscreener(dat$X, dat$Y, price = -1, verbose = FALSE), ">= 0")
  expect_error(superscreener(dat$X, dat$Y, K = 0, verbose = FALSE), "positive whole")
  expect_equal(self_sized_sets(dat$X[, "true1", drop = FALSE], dat$Y),
               setNames(rep(list("true1"), 4), SELF_SIZED))
})

test_that("user-defined learners reach the parallel workers", {
  skip_on_cran()
  SL.mine <<- function(...) SuperLearner::SL.glm(...)
  on.exit(rm(SL.mine, envir = globalenv()))
  run <- function(cores) nested_cv_superscreener(dat$X, dat$Y, prices = 0, K = 5, rankings = "lrt",
                                                 self_sized = character(0), sl_library = c("SL.mine", "SL.mean"),
                                                 V_outer = 2, seed = 1, cores = cores, verbose = FALSE)$summary
  s2 <- run(2)
  expect_equal(s2, run(1))
  expect_gt(s2$cv_auc[1], 0.6)
})
