

test_that("single realization example has local factors", {
  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  r <- 4

  result <- test_local_factors(X, r)

  expect_true(result$has_local_factors)

})


test_that("single realization example eigenvalues produced in test are same as produced outside", {

  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  M <- nrow(X)
  n <- ncol(X)
  r <- 4
  pca <- svd(X / sqrt(M))
  eig_X_outside <- pca$d^2 # eigenvalues produced in local_factors()

  eig_X_inside <- sort(eigen(t(X) %*% X / M)$values, decreasing = TRUE)

  expect_equal(eig_X_outside, eig_X_inside)

})

test_that("single realization example with eig_X missing same eig_X nonmissing", {
  skip_on_cran()
  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  M <- nrow(X)
  n <- ncol(X)
  r <- 4
  pca <- svd(X / sqrt(M))
  eig_X <- pca$d^2 # eigenvalues produced in local_factors()

  set.seed(9)
  result_inside <- test_local_factors(X, r)
  set.seed(9)
  result_outside <- test_local_factors(X, r)

  expect_equal(result_inside$Lambda, result_outside$Lambda)
  expect_equal(result_inside$has_local_factors, result_outside$has_local_factors)
  expect_equal(result_inside$n_small, result_outside$n_small)
  expect_equal(result_inside$gamma_n, result_outside$gamma_n)
  expect_equal(result_inside$h_n, result_outside$h_n)
})



test_that("single realization example with PCA matrix specified for Lambda returns no local factors", {

  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  M <- nrow(X)
  n <- ncol(X)
  r <- 4
  pca <- svd(X / sqrt(M))
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  result <- test_local_factors(X, r, Lambda0)

  expect_false(result$has_local_factors)

})

test_that("single realization with r = 3 but dim(Lambda) = (207, 4) returns error", {
  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  M <- nrow(X)
  n <- ncol(X)
  r <- 4
  pca <- svd(X / sqrt(M))
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  expect_error(test_local_factors(X, 3, Lambda0))


})







