

test_that("collating solutions does not depend on a seed, returns same result (large X)", {
  skip_on_cran()

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
  expect_equal(sol1$rotated_loadings, sol2$rotated_loadings)
  expect_equal(sol1$diagnostics$fval, sol2$diagnostics$fval)
  expect_equal(sol1$diagnostics$sol_frequency, sol2$diagnostics$sol_frequency)
})


test_that("collating solutions does not depend on seed, returns same result (single realization)", {
  skip_on_cran()

  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))
  M <- nrow(X)
  n <- ncol(X)
  r <- 4
  pca <- svd(X / sqrt(M))

  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  rotn <- find_min_rotation(Lambda0, parallel = TRUE, n_cores = 11)
  rmat_min <- rotn$R

  result1 <- collate_solutions(rotn$R, Lambda0, X)
  result2 <- collate_solutions(rotn$R, Lambda0, X)

  expect_equal(result1$rotated_loadings, result2$rotated_loadings)
  expect_equal(result1$diagnostics$R, result2$diagnostics$R)
  expect_equal(result1$diagnostics$fval, result2$diagnostics$fval)
  expect_equal(result1$diagnostics$sol_frequency, result2$diagnostics$sol_frequency)

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



test_that("nonsingular matrix columns get added to R properly", {
  r <- 4
  X <- load_matrix(testthat::test_path("fixtures", "single_realization.csv"))

  pca <- svd(X / sqrt(ncol(X)))

  Lambda0 <- sqrt(nrow(X)) * pca$v[, 1:r]

  rotn <- find_min_rotation(Lambda0, parallel = FALSE)
  rmat_min <- rotn$R

  n <- nrow(Lambda0)
  factorno <- nrow(rmat_min)
  no_randomgrid <- ncol(rmat_min)
  epsilon_rot <- 0.05

  l1_min <- colSums(abs(Lambda0 %*% rmat_min))
  sort_index <- order(l1_min)
  l1_min_sort <- l1_min[sort_index]
  rmat_min_sort <- rmat_min[, sort_index]
  rmat_min_sort <- rmat_min_sort * pracma::repmat(sign(rmat_min_sort[1, ]), factorno, 1)

  distances <- calculate_pairwise_distances(rmat_min_sort, l1_min_sort, epsilon_rot, factorno)

  candidates <- distances$rmat_min_sort %>%
    matrix_to_dataframe() %>%
    dplyr::mutate(l1_norm = distances$l1_min_sort)

  candidates <- stats::aggregate(rep(1, nrow(candidates)), by = as.list(candidates), FUN = sum) %>%
    dplyr::rename(n = x) %>%
    dplyr::arrange(l1_norm) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(non_outlier = n/gridsize(factorno) >= 0.005)

  # Construct R from the candidates
  rmat_min_unique <- candidates %>%
    dplyr::filter(non_outlier) %>%
    dplyr::select(-c(l1_norm, n, non_outlier)) %>%
    dataframe_to_matrix()

  h_n <- 1/log(n)
  amount_sparsity <- colSums(abs(Lambda0 %*% rmat_min_unique) < h_n)

  candidates <- candidates %>%
    dplyr::filter(non_outlier) %>%
    dplyr::mutate(l0_norm = amount_sparsity)

  # up until here, candidates are arranged by l1 norm
  consolidated_mins <- consolidate_local_mins(Lambda0, candidates, sorting_column = "l0_norm")

  consolidated_mins$rotated_loadings <- consolidated_mins$rotated_loadings[,1:3]

  Lambda_rotated <- fill_with_pc(consolidated_mins, Lambda0, X, r)$rotated_loadings

  is_equal_to_pca_column <- any(sapply(1:4, \(col) all(Lambda_rotated[,4] == Lambda0[,col])))

  expect_true(is_equal_to_pca_column)

})

