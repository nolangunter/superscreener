test_that("Nogueira index is 1 for identical sets and near 0 for random ones", {
  cols <- paste0("x", 1:50)
  same <- replicate(10, cols[1:5], simplify = FALSE)
  expect_equal(nogueira_index(selection_frequency(same, cols), 10), 1)
  set.seed(1)
  rand <- replicate(2000, sample(cols, 5), simplify = FALSE)
  expect_lt(abs(nogueira_index(selection_frequency(rand, cols), 2000)), 0.02)
  expect_true(is.na(nogueira_index(selection_frequency(list(cols[1]), cols), 1)))
})
