test_that("matlab and R produce same rotated Lambda in single realization, r = 2", {

  mat <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  result <- local_factors(mat, r = 2)

  result_matlab <- load_matrix(testthat::test_path("fixtures", "single_realization_lambda_rotated.csv"))
  result_matlab <- as.matrix(result_matlab)
  dimnames(result_matlab) <- NULL

  expect_equal(result$rotated_loadings, result_matlab, tolerance = 0.0001)

})

test_that("matlab and R produce same rotated Lambda in example 1 (using l1 norm to order)", {

  mat <- load_matrix(testthat::test_path("fixtures","ex1_rot_mat_matlab.csv"))
  dimnames(mat) <- NULL
  X <- load_matrix(testthat::test_path("fixtures", "example_data1.csv"))

  ex1 <- local_factors(X, 4)

  expect_equal(ex1$rotated_loadings, mat, tolerance = 0.0001)

})


test_that("local factors produces same rotation matrix (R) across same seed set outside", {
  skip_on_cran()

  mat <- load_matrix(testthat::test_path("fixtures","ex1_rot_mat_matlab.csv"))

  set.seed(916)
  result1 <- local_factors(mat, 4)

  set.seed(916)
  result2 <- local_factors(mat, 4)

  expect_identical(result1$rotation_diagnostics$R, result2$rotation_diagnostics$R)
  expect_identical(result1$rotation_diagnostics$fval, result2$rotation_diagnostics$fval)
  expect_identical(result1$rotation_diagnostics$sol_frequency, result2$rotation_diagnostics$sol_frequency)
  expect_identical(result1$rotation_diagnostics$initial_loadings, result2$rotation_diagnostics$initial_loadings)
  expect_identical(result1$rotation_diagnostics$rotated_loadings, result2$rotation_diagnostics$rotated_loadings)
  expect_identical(result1$rotation_diagnostics$initial_loadings, result2$rotation_diagnostics$initial_loadings)

})



test_that("matlab and R produce same rotated Lambda in example 1 (using l0 norm to order)", {
  skip_on_cran()

  mat <- load_matrix(testthat::test_path("fixtures","lambda_rotated_ex1.csv"))
  dimnames(mat) <- NULL
  X <- load_matrix(testthat::test_path("fixtures", "example_data1.csv"))

  ex1 <- local_factors(X, 4)

  expect_equal(ex1$rotated_loadings, mat, tolerance = 0.01)

})

test_that("matlab and csv R produce same rotated Lambda in example 1 (using l0 norm to order)", {
  skip_on_cran()

  mat <- load_matrix(testthat::test_path("fixtures","lambda_rotated_ex1.csv"))
  dimnames(mat) <- NULL
  X <- read.csv(testthat::test_path("fixtures", "example_data1.csv"), header = FALSE)

  ex1 <- local_factors(X, 4)

  expect_equal(ex1$rotated_loadings, mat, tolerance = 0.01)

})

test_that("local_factors returns error when missing values are in the data matrix", {

  X <- load_matrix(testthat::test_path("fixtures", "example_data1.csv"))
  indices <- sample(1:length(X), 100)
  X[indices] <- NA
  expect_error(local_factors(X, 4))

})



