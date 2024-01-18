find_local_factors <- function(X, r, Lambda0) {
  # Function to find the sparsest rotation of the leading Principal Components
  # of a T by n matrix X. 
  # Under sparsity in the loading matrix this will identify the true loading matrix.
  #
  # returns two arguments:
  #   Lambda0: If not provided, Principal Component estimate
  #   Lambda_rotated: Rotation of loading matrix with smallest l1-norm
  T <- nrow(X)
  n <- ncol(X)
  svd_X <- svd(X / sqrt(T))
  eig_X <- svd_X$d^2
  if (missing(Lambda0)) {
    Lambda0 <- sqrt(n) * svd_X$v[, 1:r]
  }
  # compute the rotated solution with minimal l1-norm
  rmat_min_results <- find_min_rotation(Lambda0) #Finds solution for each point in grid
  rmat_min <- rmat_min_results$R
  rotation_results <- collate_solutions(rmat_min, Lambda0, eig_X) #Combine into candidates
  Lambda_rotated <- rotation_results$Lambda_rotated
  diagnostics <- rotation_results$diagnostics
  return(list(Lambda0 = Lambda0, Lambda_rotated = Lambda_rotated, diagnostics = diagnostics))
}


