
test_that("find_local_factors() returns same result with same seed with and without Lambda0 argument, small example", {
  X <- readr::read_csv(testthat::test_path("fixtures", "single_realization.csv"), col_names = FALSE)
  X <- as.matrix(X)
  r <- 2
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  eig_X <- pca$d^2
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result <- find_local_factors(X, r)
  set.seed(916)
  with_lambda_result <- find_local_factors(X, r, Lambda0)

  expect_equal(result$diagnostics$R, with_lambda_result$diagnostics$R, tolerance = 0.0001)

})


test_that("find_local_factors() returns same result with same seed with and without Lambda0 argument, larger example", {
  X <- readr::read_csv(testthat::test_path("fixtures", "example_data1.csv"), col_names = FALSE) %>%
    as.matrix()
  X <- as.matrix(X)
  r <- 4
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  eig_X <- pca$d^2
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result <- find_local_factors(X, r)
  set.seed(916)
  with_lambda_result <- find_local_factors(X, r, Lambda0)

  expect_equal(result$diagnostics$R, with_lambda_result$diagnostics$R, tolerance = 0.0001)

})




