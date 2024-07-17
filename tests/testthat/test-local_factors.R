test_that("matlab and R produce same rotated Lambda in single realization, r = 2", {

  mat <- readr::read_csv(testthat::test_path("fixtures", "single_realization.csv"), col_names = FALSE)
  mat <- as.matrix(mat)
  result <- local_factors(mat, r = 2)

  result_matlab <- readr::read_csv(testthat::test_path("fixtures", "single_realization_lambda_rotated.csv"), col_names = FALSE)
  result_matlab <- as.matrix(result_matlab)
  dimnames(result_matlab) <- NULL

  expect_equal(result$Lambda_rotated, result_matlab, tolerance = 0.0001)

})

test_that("matlab and R produce same rotated Lambda in example 1 (using l1 norm to order)", {

  mat <- readr::read_csv(testthat::test_path("fixtures","ex1_rot_mat_matlab.csv"), col_names = FALSE) %>%
    as.matrix()
  dimnames(mat) <- NULL
  X <- readr::read_csv(testthat::test_path("fixtures", "example_data1.csv"), col_names = FALSE) %>%
    as.matrix()

  ex1 <- local_factors(X, 4)

  expect_equal(ex1$Lambda_rotated, mat, tolerance = 0.0001)

})


test_that("local factors produces same rotation matrix (R) across same seed set outside", {
  mat <- readr::read_csv(testthat::test_path("fixtures","ex1_rot_mat_matlab.csv"), col_names = FALSE) %>%
    as.matrix()

  set.seed(916)
  result1 <- local_factors(mat, 4)

  set.seed(916)
  result2 <- local_factors(mat, 4)

  expect_equal(result1$rotation_diagnostics$R, result2$rotation_diagnostics$R, tolerance = 0.0001)

})



test_that("matlab and R produce same rotated Lambda in example 1 (using l0 norm to order)", {

  mat <- readr::read_csv(testthat::test_path("fixtures","lambda_rotated_ex1.csv"), col_names = FALSE) %>%
    as.matrix()
  dimnames(mat) <- NULL
  X <- readr::read_csv(testthat::test_path("fixtures", "example_data1.csv"), col_names = FALSE) %>%
    as.matrix()

  ex1 <- local_factors(X, 4)

  expect_equal(ex1$Lambda_rotated, mat, tolerance = 0.01)

})



