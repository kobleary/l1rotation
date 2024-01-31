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


test_that("negative input to spherical_to_cartesian returns error", {

})







