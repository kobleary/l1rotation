
test_that("find_local_factors() returns same result with same seed with and without Lambda0 argument, small example", {
  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  r <- 2
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result <- find_local_factors(X, r)
  set.seed(916)
  with_lambda_result <- find_local_factors(X, r, Lambda0)

  expect_equal(result$diagnostics$R, with_lambda_result$diagnostics$R)
  expect_equal(result$diagnostics$fval, with_lambda_result$diagnostics$fval)
  expect_equal(result$diagnostics$sol_frequency, with_lambda_result$diagnostics$sol_frequency)
  expect_equal(result$Lambda0, with_lambda_result$Lambda0)
  expect_equal(result$Lambda, with_lambda_result$Lambda)
})


test_that("find_local_factors() returns same result with same seed with and without Lambda0 argument, larger example", {
  skip_on_cran()

  X <- load_matrix(testthat::test_path("fixtures", "example_data1.csv"))
  r <- 4
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result <- find_local_factors(X, r)
  set.seed(916)
  with_lambda_result <- find_local_factors(X, r, Lambda0)

  expect_equal(result$diagnostics$R, with_lambda_result$diagnostics$R)
  expect_equal(result$diagnostics$fval, with_lambda_result$diagnostics$fval)
  expect_equal(result$diagnostics$sol_frequency, with_lambda_result$diagnostics$sol_frequency)
  expect_equal(result$Lambda0, with_lambda_result$Lambda0)
  expect_equal(result$Lambda, with_lambda_result$Lambda)

})

test_that("find_local_factors() returns same result with same seed with and without Lambda0 argument, randomly generated data", {
  skip_on_cran()

  X <- matrix(stats::rnorm(100*78), nrow = 100, ncol = 78)

  r <- 4
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result <- find_local_factors(X, r)
  set.seed(916)
  with_lambda_result <- find_local_factors(X, r, Lambda0)

  expect_equal(result$diagnostics$R, with_lambda_result$diagnostics$R)
  expect_equal(result$diagnostics$fval, with_lambda_result$diagnostics$fval)
  expect_equal(result$diagnostics$sol_frequency, with_lambda_result$diagnostics$sol_frequency)
  expect_equal(result$Lambda0, with_lambda_result$Lambda0)
  expect_equal(result$Lambda, with_lambda_result$Lambda)


})




