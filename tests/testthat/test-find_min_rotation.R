# Should not include library() calls

test_that("sphere_to_cart o cart_to_sphere = matrix, 2 factors (2 x 3), ", {

  mat <- normalize(
    matrix(
     c(-0.3896802, 0.1362248, -0.2546876,
      0.9209502, -0.9906779,  0.9670234),
    nrow = 2)
    )

  expect_equal(spherical_to_cartesian(cartesian_to_spherical(mat)), mat)

})

test_that("cart_to_sphere o sphere_to_cart = matrix, 2 factors (2 x 3), ", {

  mat <- normalize(
    matrix(
      c(0.3896802, 0.1362248,  0.2546876,
        0.9209502, 0.9906779,  0.9670234),
      nrow = 2)
  )

  expect_equal(cartesian_to_spherical(spherical_to_cartesian(mat)), mat)

})



test_that("sphere_to_cart o cart_to_sphere = matrix, 3 factors (3 x 3), ", {

  mat <- normalize(
    matrix(
      c(0.1232947, 1.6504757, 0.1633116,
      0.8283809, 0.3916161, 0.3547165,
      0.3054703, 0.5366281, 0.9906242),
      nrow = 3)
  )

  expect_equal(spherical_to_cartesian(cartesian_to_spherical(mat)), mat)

})

test_that("cart_to_sphere o sphere_to_cart = matrix, 3 factors (3 x 3), ", {

  mat <- normalize(
    matrix(
      c(0.1232947, 1.6504757, 0.1633116,
        0.8283809, 0.3916161, 0.3547165,
        0.3054703, 0.5366281, 0.9906242),
      nrow = 3)
  )

  expect_equal(cartesian_to_spherical(spherical_to_cartesian(mat)), mat)

})


test_that("normalization works", {

  norm_mat <- normalize(matrix(stats::rnorm(3 * 4), nrow = 3), p = 2)

  # Calculate norms by hand
  norms <- purrr::map_dbl(1:nrow(norm_mat), ~ sqrt(sum(norm_mat[, .]^2)))

  expect_equal(norms, rep(1, nrow(norm_mat)))

})


test_that("matlab and R produce same rotated Lambda in example 1", {

  mat <- readr::read_csv(testthat::test_path("fixtures","ex1_rot_mat_matlab.csv"), col_names = FALSE) %>%
    as.matrix()
  dimnames(mat) <- NULL
  X <- readr::read_csv(testthat::test_path("fixtures", "example_data1.csv"), col_names = FALSE) %>%
    as.matrix()

  ex1 <- local_factors(X, 4)

  expect_equal(ex1$Lambda_rotated, mat, tolerance = 0.0001)

})

test_that("objective function works with spherical_to_cartesian(), r = 3", {

  r <- 3
  Lambda <- matrix(stats::rnorm(r * 100), nrow = 100)
  theta <- matrix(stats::rnorm(r * 4), nrow = r) %>% cartesian_to_spherical()
  theta <- theta[,1]

  computed <- objectivefcn_spherical(theta, Lambda)

  R <- rep(0, r)
  R[1] <- cos(theta[1])
  if(r > 2){
    for (kk in 2:(r-1)) {
      R[kk] <- prod(sin(theta[1:(kk-1)])) * cos(theta[kk])
    }
  }
  R[r] <- prod(sin(theta))

  expected <- sum(abs(Lambda %*% R))

  expect_equal(computed, expected, tolerance = 0.001)


})

test_that("objective function works with spherical_to_cartesian(), r = 2", {

  r <- 2
  Lambda <- matrix(stats::rnorm(r * 100), nrow = 100)
  theta <- matrix(stats::rnorm(r * 4), nrow = r) %>% cartesian_to_spherical()
  theta <- theta[,1]

  computed <- objectivefcn_spherical(theta, Lambda)

  R <- rep(0, r)
  R[1] <- cos(theta[1])
  if(r > 2){
    for (kk in 2:(r-1)) {
      R[kk] <- prod(sin(theta[1:(kk-1)])) * cos(theta[kk])
    }
  }
  R[r] <- prod(sin(theta))

  expected <- sum(abs(Lambda %*% R))

  expect_equal(computed, expected, tolerance = 0.001)

})

test_that("objective function with length(theta) =/= ncol(Lambda) returns non-conformable error", {

  r <- 2
  Lambda <- matrix(stats::rnorm((r + 1) * 100), nrow = 100)
  theta <- matrix(stats::rnorm((r-1) * 4), nrow = r) %>% cartesian_to_spherical()
  theta <- theta[,1]

  expect_error(objectivefcn_spherical(theta, Lambda), "non-conformable")

})


test_that("negative input to spherical_to_cartesian returns error", {

})







