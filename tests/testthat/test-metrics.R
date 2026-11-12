test_that("auc matches the pairwise definition", {
  set.seed(1); y <- rbinom(200, 1, 0.3); p <- runif(200)
  pairs <- outer(p[y == 1], p[y == 0], ">") + 0.5 * outer(p[y == 1], p[y == 0], "==")
  expect_equal(auc(y, p), mean(pairs))
})

test_that("cv_auc averages fold AUCs and a paired difference with itself is zero", {
  set.seed(2); y <- rbinom(300, 1, 0.4); p <- plogis(rnorm(300) + y); f <- stratified_folds(y, 5)
  a <- cv_auc(y, p, f)
  expect_equal(a$estimate, mean(sapply(1:5, function(v) auc(y[f == v], p[f == v]))))
  expect_true(a$se > 0 && a$lower < a$estimate && a$upper > a$estimate)
  d <- cv_auc(y, p, f, p_ref = p)
  expect_equal(d$estimate, 0); expect_equal(d$se, 0)
})

test_that("brier and neg_loglik", {
  expect_equal(brier(c(0, 1), c(0.2, 0.6)), (0.04 + 0.16) / 2)
  expect_true(is.finite(neg_loglik(c(0, 1), c(0, 1))))
})

test_that("stratified folds balance the outcome", {
  set.seed(3); y <- rbinom(500, 1, 0.3); f <- stratified_folds(y, 5)
  expect_equal(sort(unique(f)), 1:5)
  expect_true(max(table(f, y)[, "1"]) - min(table(f, y)[, "1"]) <= 1)
})
