#' Find the most sparse rotation of the leading Principal Components of a t by n matrix X.
#'
#' @inheritParams local_factors
#' @param Lambda0 A matrix that represents an orthonormal basis of the loading space. If not supplied, PCA is used.
#'
#' @returns Returns a list with the following components:
#'  * `Lambda0` (if not supplied) the principal component estimate of the loading matrix (orthonormal)
#'  * `Lambda_rotated` a matrix that is the rotation of the loading matrix that produces the smallest l1-norm.
#'  * `rotation_diagnostics` a list containing 3 components"
#'      * `R` the rotation matrix that when used to rotate `Lambda0` produces the smallest l1-norm.
#'      * `l1_norm` a vector of length `r` containing the value of the l1 norm each solution generates
#'      * `sol_frequency` a vector of length `r` containing the frequency in the initial grid of each solution
#'
#' @export
find_local_factors <- function(X, r, Lambda0) {
  # Function to find the sparsest rotation of the leading Principal Components
  # of a M by n matrix X.
  # Under sparsity in the loading matrix this will identify the true loading matrix.
  #
  # returns three arguments:
  #   Lambda0: Principal Component estimate
  #   Lambda_rotated: Rotation of loading matrix with smallest l1-norm
  #   diagnostics: A list of diagnostics

  # X cannot have missing values
  stopifnot(!any(is.na(X)))
  stopifnot(!any(is.infinite(X)))

  M <- nrow(X)
  n <- ncol(X)
  svd_X <- svd(X / sqrt(M), nu = M, nv = n)
  eig_X <- svd_X$d^2
  if (missing(Lambda0)) {
    Lambda0 <- sqrt(n) * svd_X$v[, 1:r]
  }
  # compute the rotated solution with minimal l1-norm
  rmat_min_results <- find_min_rotation(Lambda0) #Finds solutions across a large grid of starting points
  rmat_min <- rmat_min_results$R
  rotation_results <- collate_solutions(rmat_min, Lambda0, X) #Combine large number of solutions into candidates
  Lambda_rotated <- rotation_results$Lambda_rotated
  diagnostics <- rotation_results$diagnostics
  return(list(Lambda0 = Lambda0, Lambda_rotated = Lambda_rotated, diagnostics = diagnostics))
}


