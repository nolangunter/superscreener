d <- make_design(p = 40, mc_n = 2e4)
dat <- simulate_data(200, d, seed = 1)

test_that("design has the documented roles and calibration", {
  expect_equal(as.vector(table(d$roles)[c("true", "proxy", "noise")]), c(10, 20, 10))
  expect_equal(d$binary$prev, 0.3)
  expect_equal(ncol(dat$X), 40); expect_true(all(dat$Y %in% 0:1))
})

test_that("rankings are permutations, constant columns go last", {
  X <- dat$X; X$const <- 1
  rk <- rank_columns(X, dat$Y)
  expect_equal(colnames(rk), RANKINGS)
  expect_equal(unname(sort(rk[, "consensus"][-41])), 1:40)
  expect_true(all(rk["const", c("lrt", "lasso", "rf", "consensus")] == 41))
  expect_equal(top_k(rk[, "lrt"], 3), names(sort(rk[, "lrt"]))[1:3])
})

test_that("self-sized screeners return column names", {
  s <- self_sized_sets(dat$X, dat$Y, seed = 1)
  expect_named(s, SELF_SIZED)
  expect_true(all(unlist(s) %in% colnames(dat$X)))
  expect_gte(length(s$sl_glmnet), 2)
})

test_that("expand_factors charges a categorical variable once", {
  df <- data.frame(age = c(1, 2, 3, 4), grp = c("a", "b", "c", "a"), flag = c(TRUE, FALSE, TRUE, TRUE))
  ex <- expand_factors(df)
  expect_equal(unname(ex$groups), c("age", "grp", "grp", "flag"))
  expect_true(all(vapply(ex$X, is.numeric, TRUE)))
  expect_error(expand_factors(data.frame(a = c(1, NA))), "missing")
})
