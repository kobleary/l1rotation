

test_that("collating solutions does not depend on a seed, returns same result (large X)", {

  initial_draws <- load_matrix(testthat::test_path("fixtures", "initial_draws_ex1.csv"))
  X <- load_matrix(testthat::test_path("fixtures", "example_data1.csv"))

  Lambda <- load_matrix(testthat::test_path("fixtures", "Lambda_ex1.csv"))

  r <- 4
  M <- nrow(X)
  n <- ncol(X)

  result <- find_min_rotation(Lambda)
  rmat_min <- result$R

  sol1 <- collate_solutions(rmat_min, Lambda, X)
  sol2 <- collate_solutions(rmat_min, Lambda, X)

  expect_equal(sol1$diagnostics$R, sol2$diagnostics$R)

})


test_that("collating solutions does not depend on seed, returns same result (single realization)", {
  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  M <- nrow(X)
  n <- ncol(X)
  r <- 4
  pca <- svd(X / sqrt(M))

  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  rotn <- find_min_rotation(Lambda0, parallel = TRUE)
  rmat_min <- rotn$R

  result1 <- collate_solutions(rotn$R, Lambda0, X)
  result2 <- collate_solutions(rotn$R, Lambda0, X)

  expect_equal(result1$R, result2$R)

})


test_that("||a - b||^2 is equal to ||a||^2 + ||b||^2 - 2<a, b> for an entry from a random matrix", {

  rmat_min_sort <- matrix(stats::rnorm(200 * 4), ncol = 200)

  cross_prod <- crossprod(rmat_min_sort)
  diag_vals <- diag(cross_prod)
  ind <- sample(1:ncol(rmat_min_sort), 2)

  norm_slow <- pracma::Norm(rmat_min_sort[,ind[1]] - rmat_min_sort[, ind[2]], p = 2)
  norm_matrix <- sqrt(diag_vals[ind[1]] + diag_vals[ind[2]] - 2 * cross_prod[ind[1],ind[2]])

  expect_equal(norm_slow, norm_matrix)
})


