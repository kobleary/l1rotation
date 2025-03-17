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
  norms <- sapply(1:nrow(norm_mat), \(row) sqrt(sum(norm_mat[, row]^2)))

  expect_equal(norms, rep(1, nrow(norm_mat)))

})


test_that("data.frame to find_min_rotation returns an error", {
  r <- 4
  Lambda <- matrix(stats::rnorm(r * 100), nrow = 100) |>
    as.data.frame()
  expect_error(find_min_rotation(Lambda))

})

test_that("missing values in Lambda returns an error", {
  r <- 8
  Lambda <- matrix(stats::rnorm(r * 100), nrow = 100)

  na_indexes <- sample(1:length(Lambda), 10)
  Lambda[na_indexes] <- NA

  expect_error(find_min_rotation(Lambda, parallel = FALSE), "Lambda contains missing or infinite values.")

})

test_that("non-numeric values in Lambda returns an error", {
  r <- 8
  Lambda <- matrix(stats::rnorm(r * 100), nrow = 100)

  indexes <- sample(1:length(Lambda), 10)
  Lambda[indexes] <- "string"

  expect_error(find_min_rotation(Lambda, parallel = FALSE), "Lambda contains non-numeric values.")
})



test_that("objective function works with spherical_to_cartesian(), r = 3", {

  r <- 3
  Lambda <- matrix(stats::rnorm(r * 100), nrow = 100)
  theta <- matrix(stats::rnorm(r * 4), nrow = r) |> cartesian_to_spherical()
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

  expect_equal(computed, expected)


})

test_that("objective function works with spherical_to_cartesian(), r = 2", {

  r <- 2
  Lambda <- matrix(stats::rnorm(r * 100), nrow = 100)
  theta <- matrix(stats::rnorm(r * 4), nrow = r) |> cartesian_to_spherical()
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

  expect_equal(computed, expected)

})

test_that("objective function with length(theta) =/= ncol(Lambda) returns non-conformable error", {

  r <- 2
  Lambda <- matrix(stats::rnorm((r + 1) * 100), nrow = 100)
  theta <- matrix(stats::rnorm((r-1) * 4), nrow = r) |> cartesian_to_spherical()
  theta <- theta[,1]

  expect_error(objectivefcn_spherical(theta, Lambda), "non-conformable")

})



test_that("single realization returns same R matrix", {
  initial_draws <- load_matrix(testthat::test_path("fixtures", "initial_draws_ex1.csv"))
  X <- load_matrix(testthat::test_path("fixtures", "example_data1.csv"))
  Lambda <- load_matrix(testthat::test_path("fixtures", "Lambda_ex1.csv"))

  r <- 4
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M))

  r <- ncol(Lambda)
  no_draws <- ncol(initial_draws)
  l1_norm <- rep(0, no_draws)
  exitflag <- rep(0, no_draws)

  theta <- cartesian_to_spherical(initial_draws)

  # Optimization in polar coordinates happens w.r.t. theta
  angles <- theta
  l <- nrow(angles)
  for (rep in 1:no_draws) {
    starting_point <- theta[, rep]
    result <- stats::optim(
      starting_point,
      objectivefcn_spherical, Lambda = Lambda,
      control = list(maxit = 200 * l, ndeps = 1e-8, reltol = 1e-8, warn.1d.NelderMead = FALSE),
      method = 'Nelder-Mead'
    )

    angles[, rep] <- result$par
    l1_norm[rep] <- result$value
    exitflag[rep] <- result$convergence
  }

  # Convert back to cartesian coordinates, need to edit to generalize across minimum number of factors (requires at least 2)
  R <- spherical_to_cartesian(angles)

  result <- collate_solutions(R, Lambda, X)


  R_matlab <- load_matrix(testthat::test_path("fixtures", "R_ex1.csv"))
  result_matlab <- collate_solutions(R_matlab, Lambda, X)

  expect_equal(result$diagnostics$R, result_matlab$diagnostics$R, tolerance = 0.0001)
  expect_equal(result$Lambda_rotated, result_matlab$Lambda_rotated, tolerance = 0.0001)
  expect_equal(result$fval, result_matlab$fval, tolerance = 0.0001)
  expect_equal(result$sol_frequency, result_matlab$sol_frequency, tolerance = 0.0001)

})

test_that("find_min_rotation(), check that parallel gives same R with same seed (4 factors)", {
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
  result_1 <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  set.seed(916)
  result_2 <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)

  expect_equal(result_1$R, result_2$R)
  expect_equal(result_1$exitflag, result_2$exitflag)
  expect_equal(result_1$l1_norm, result_2$l1_norm)


})



test_that("find_min_rotation(), check that non-parallel gives same R with same seed (4 factors)", {
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
  result_1 <- find_min_rotation(Lambda0, parallel = FALSE)
  set.seed(916)
  result_2 <- find_min_rotation(Lambda0, parallel = FALSE)

  expect_equal(result_1$R, result_2$R)
  expect_equal(result_1$exitflag, result_2$exitflag)
  expect_equal(result_1$l1_norm, result_2$l1_norm)

})


test_that("find_min_rotation(), check that non-parallel/parallel gives same R with same seed (4 factors)", {
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
  result_par <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  set.seed(916)
  result_nonpar <- find_min_rotation(Lambda0, parallel = FALSE)

  expect_equal(result_par$R, result_nonpar$R)
  expect_equal(result_par$exitflag, result_nonpar$exitflag)
  expect_equal(result_par$l1_norm, result_nonpar$l1_norm)

})

test_that("find_min_rotation(), check that non-parallel/parallel gives same R with same seed (8 factors)", {
  skip_on_cran()

  X <- load_matrix(testthat::test_path("fixtures", "example_data2.csv"))
  r <- 8
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result_par <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  set.seed(916)
  result_nonpar <- find_min_rotation(Lambda0, parallel = FALSE)

  expect_equal(result_par$R, result_nonpar$R)
  expect_equal(result_par$exitflag, result_nonpar$exitflag)
  expect_equal(result_par$l1_norm, result_nonpar$l1_norm)
})



test_that("find_min_rotation(), check that parallel gives same R with same seed (4 factors)", {
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
  result_1 <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  set.seed(916)
  result_2 <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)

  expect_equal(result_1$R, result_2$R)
  expect_equal(result_1$exitflag, result_2$exitflag)
  expect_equal(result_1$l1_norm, result_2$l1_norm)


})



test_that("find_min_rotation(), check that non-parallel gives same R with same seed (4 factors)", {
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
  result_1 <- find_min_rotation(Lambda0, parallel = FALSE)
  set.seed(916)
  result_2 <- find_min_rotation(Lambda0, parallel = FALSE)

  expect_equal(result_1$R, result_2$R)
  expect_equal(result_1$exitflag, result_2$exitflag)
  expect_equal(result_1$l1_norm, result_2$l1_norm)

})


test_that("find_min_rotation(), check that non-parallel/parallel gives same R with same seed (4 factors)", {
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
  result_par <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  set.seed(916)
  result_nonpar <- find_min_rotation(Lambda0, parallel = FALSE)

  expect_equal(result_par$R, result_nonpar$R)
  expect_equal(result_par$exitflag, result_nonpar$exitflag)
  expect_equal(result_par$l1_norm, result_nonpar$l1_norm)

})

# Smaller examples with single_realization ------------------------------------

test_that("find_min_rotation(), check that parallel gives same R with same seed (small example)", {

  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  r <- 2
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result_1 <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  set.seed(916)
  result_2 <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)

  expect_equal(result_1$R, result_2$R)
  expect_equal(result_1$exitflag, result_2$exitflag)
  expect_equal(result_1$l1_norm, result_2$l1_norm)

})



test_that("find_min_rotation(), check that non-parallel gives same R with same seed (small example)", {

  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  r <- 2
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result_1 <- find_min_rotation(Lambda0, parallel = FALSE)
  set.seed(916)
  result_2 <- find_min_rotation(Lambda0, parallel = FALSE)

  expect_equal(result_1$R, result_2$R)
  expect_equal(result_1$exitflag, result_2$exitflag)
  expect_equal(result_1$l1_norm, result_2$l1_norm)


})


test_that("find_min_rotation(), check that non-parallel/parallel gives same R with same seed (4 factors)", {

  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  r <- 2
  M <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(M), nu = M, nv = n)
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Find minimum rotation, test for local factors
  set.seed(916)
  result_par <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  set.seed(916)
  result_nonpar <- find_min_rotation(Lambda0, parallel = FALSE)

  expect_equal(result_par$R, result_nonpar$R)
  expect_equal(result_par$exitflag, result_nonpar$exitflag)
  expect_equal(result_par$l1_norm, result_nonpar$l1_norm)


})


