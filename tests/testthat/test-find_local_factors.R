
max_cosine_similarity <- function(est_mat, true_mat){
  colMaxs(abs(t(est_mat) %*% true_mat)/nrow(true_mat))
}

orthonormalize <- function(matrix){
  r <- ncol(matrix)
  n <- nrow(matrix)


  orthonorm_matrix <- sqrt(n) * pracma::gramSchmidt(matrix)$Q

  return(orthonorm_matrix)

}

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


test_that("find_local_factors() returns matrix similar to truth using a orthonormal rotation of the true Lambda", {

  i <- 10

  X <- load_matrix(testthat::test_path("fixtures", stringr::str_glue("X_{i}.csv")))
  true_lambda <- load_matrix(testthat::test_path("fixtures", stringr::str_glue("true_lambda_{i}.csv")))
  e <- matrix(rnorm(n * r, sd = 0.1), nrow = n, ncol = r)

  true_lambda_rotated_plus_noise <- orthonormalize(true_lambda + e)

  M <- nrow(X)
  n <- ncol(X)
  r <- 4

  lambda_rotated_true_plus_noise <- find_local_factors(X, r, true_lambda_rotated_plus_noise)$Lambda_rotated
  cos_sim_true_plus_noise <- max_cosine_similarity(lambda_rotated_true_plus_noise, true_lambda)[-1]
  expect_true(all(cos_sim_true_plus_noise > 0.98))

})





