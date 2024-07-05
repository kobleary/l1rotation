collate_solutions <- function(rmat_min, Lambda_0, eig_X) {

  tictoc::tic("collate_solutions")

  n <- nrow(Lambda_0)
  factorno <- nrow(rmat_min)
  no_randomgrid <- ncol(rmat_min)
  epsilon_rot <- 0.05

  l1_min <- colSums(abs(Lambda_0 %*% rmat_min))
  sort_index <- order(l1_min)
  l1_min_sort <- l1_min[sort_index]
  rmat_min_sort <- rmat_min[, sort_index]
  rmat_min_sort <- rmat_min_sort * pracma::repmat(sign(rmat_min_sort[1, ]), factorno, 1)

  tictoc::tic() # This step takes 44 seconds with 8 factors and 3000 draws (probably can be improved)
  norms <- matrix(0, nrow = no_randomgrid, ncol = no_randomgrid)
  for (i in 1:no_randomgrid) {
    for (j in i:no_randomgrid) {

      norms[i,j] <- (pracma::Norm(rmat_min_sort[, i] - rmat_min_sort[, j], p = 2) / sqrt(factorno))
      if (norms[i,j] < epsilon_rot) {
        rmat_min_sort[, j] <- rmat_min_sort[, i]
        l1_min_sort[j] <- l1_min_sort[i]
      }
    }
  }

  candidates <- rmat_min_sort %>%
    matrix_to_dataframe() %>%
    dplyr::mutate(l1_norm = l1_min_sort) %>%
    dplyr::group_by(dplyr::across(tidyselect::everything())) %>%
    dplyr::count() %>%
    dplyr::arrange(l1_norm) %>%
    dplyr::ungroup() %>%
    mutate(non_outlier = n/gridsize(factorno) >= 0.005)


  # Construct R from the candidates
  rmat_min_unique <- candidates %>%
    filter(non_outlier) %>%
    dplyr::select(-c(l1_norm, n, non_outlier)) %>%
    dataframe_to_matrix()

  h_n <- 1/log(n)
  amount_sparsity <- colSums(abs(Lambda_0 %*% rmat_min_unique) < h_n)

  candidates <- candidates %>%
    filter(non_outlier) %>%
    mutate(l0_norm = amount_sparsity)

  # up until here, candidates are arranged by l1 norm
  consolidated_mins <- consolidate_local_mins(Lambda_0, candidates, sorting_column = "l0_norm")

  Lambda_rotated <- consolidated_mins$Lambda_rotated
  R <- consolidated_mins$R

  # Fill with PCs far from collinear if number of candidates < number of factors
  candidateno <- ncol(Lambda_rotated)
  I <- diag(factorno)
  l1_norm_update <- c()

  while (candidateno < factorno) {
    consolidated_mins <- fill_with_pc(consolidated_mins, factorno, I, Lambda_0, eig_X, R, l1_norm_update)
    candidateno <- ncol(consolidated_mins$Lambda_rotated)
    l1_norm_update <- consolidated_mins$l1_norm_update
    R <- consolidated_mins$R
  }

  loadings <- matrix_to_dataframe(R) %>%
    dplyr::left_join(candidates) %>%
    dplyr::mutate(
      l1_norm = c(l1_norm[!is.na(l1_norm)], l1_norm_update),
      n = ifelse(is.na(n), 0, n)
    )

  tictoc::toc()

  return(
    list(
      diagnostics = list(R = R, fval = loadings$l1_norm, sol_frequency = loadings$n),
      Lambda_rotated = consolidated_mins$Lambda_rotated
    )
  )

}

dataframe_to_matrix <- function(dataframe) {
  dataframe %>% as.matrix() %>% t()
}

matrix_to_dataframe <- function(matrix){
  matrix %>% t() %>% as.data.frame()
}

consolidate_local_mins <- function(Lambda_0, candidates, sorting_column = "l0_norm") {

  # Sorting column determines what is used to pick local minima
  # l1_norm
  # n - number of occurences
  # l0_norm - approximate l0-norm

  if (sorting_column != "l1_norm") {

    candidates <- candidates %>%
      dplyr::arrange(dplyr::across(tidyselect::all_of(sorting_column), desc))

    rmat_min_unique <- candidates %>%
      dplyr::select(tidyr::starts_with("V")) %>%
      dataframe_to_matrix()
  }

  # Consolidate candidates
  Lambda_rotated <- Lambda_0 %*% rmat_min_unique[, 1] # First candidate
  R <- rmat_min_unique[, 1]

  factorno <- nrow(rmat_min_unique)
  upperK <- ncol(rmat_min_unique)
  n <- length(Lambda_0)

   for (kk in 2:upperK) {

     temp <- cbind(Lambda_rotated, Lambda_0 %*% rmat_min_unique[, kk])
     not_singular <- min(eigen(t(temp) %*% temp)$values) / n > sqrt(1 / factorno) / 3 &&
       min(eigen(t(temp) %*% temp)$values) / n > min(eigen(t(Lambda_rotated) %*% Lambda_rotated)$values) / n / 4

     if(not_singular){
       Lambda_rotated <- temp
       R <- cbind(R, rmat_min_unique[, kk])
     }
   }

  return(list(Lambda_rotated = Lambda_rotated, R = R))

}

fill_with_pc <- function(consolidated_mins, factorno, I, Lambda_0, eig_X, R, l1_norm_update){
  Lambda_rotated <- consolidated_mins$Lambda_rotated

  min_eig <- rep(0, factorno)
  for (ell in 1:factorno) {
    temp <- cbind(Lambda_rotated, Lambda_0[, ell])
    min_eig[ell] <- min(eigen(t(temp) %*% temp)$values)
  }
  index <- which.min(min_eig * sqrt(eig_X[1:factorno]))
  print(index)
  print(min_eig)
  Lambda_rotated <- cbind(Lambda_rotated, Lambda_0[, index])
  R <- cbind(R, I[, index])
  l1_norm_update <- c(l1_norm_update, sum(abs(Lambda_0[, index])))

  return(list(Lambda_rotated = Lambda_rotated, R = R, l1_norm_update = l1_norm_update))

}
