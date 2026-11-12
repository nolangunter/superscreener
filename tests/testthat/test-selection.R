cand <- data.frame(candidate = c("a_K5", "a_K10", "b", "c_K20"), K = c(5, 10, NA, 20),
                   size = c(5, 10, 3, 20), cv_risk = c(0.170, 0.165, 0.180, 0.164))

test_that("price 0 is the discrete Super Learner", {
  expect_equal(choose_candidate(cand, 0, 0.21)$candidate, "c_K20")
})

test_that("a higher price picks a smaller set, and lambda follows the formula", {
  ch <- choose_candidate(cand, 0.01, 0.21)
  expect_equal(ch$lambda, 0.01 * (0.21 - 0.164))
  expect_equal(ch$candidate, "a_K10")
  expect_equal(choose_candidate(cand, 0.2, 0.21)$candidate, "b")
})

test_that("no price when nothing beats the intercept; ties go to the smaller set", {
  expect_equal(choose_candidate(cand, 0.5, 0.1)$lambda, 0)
  tie <- data.frame(candidate = c("x_K10", "y"), K = c(10, NA), size = c(10, 4), cv_risk = c(0.2, 0.2))
  expect_equal(choose_candidate(tie, 0, 0.3)$candidate, "y")
})

test_that("min_size drops small candidates", {
  expect_equal(choose_candidate(cand, 0.1, 0.21, min_size = 5)$candidate, "a_K5")
})
