#' Find the rotation of the loading matrix with the smallest l1-norm, as in [local_factors()], with additional flexibility.
#'
#'
#' @description
#' Find the most sparse rotation of an orthonormal basis of the loading space of a t by n matrix X. Additional flexibility with the `initial_loadings` argument allows the user to specify any orthonormal basis rather than defaulting to PCA.
#'
#' @inheritParams local_factors
#' @param initial_loadings Matrix that represents an orthonormal basis of the loading space. If not supplied, PCA is used by default in this function and also in `local_factors`.
#'
#' @returns Returns a list with the following components:
#'  * `initial_loadings` Principal Component estimate of the loading matrix (if not supplied).
#'  * `rotated_loadings` Matrix that is the rotation of the loading matrix that produces the smallest l1-norm.
#'  * `rotation_diagnostics` A list containing 3 components:
#'      * `R` Rotation matrix that when used to rotate `initial_loadings` produces the smallest l1-norm.
#'      * `l1_norm` Vector of length `r` containing the value of the l1 norm each solution generates.
#'      * `sol_frequency` Vector of length `r` containing the frequency in the initial grid of each solution.
#'
#' @export
#'
#' @examples
#' # Minimal example with 2 factors, where X is a 224 by 207 matrix
#' r <- 2
#' M <- nrow(example_data)
#' n <- ncol(example_data)
#'
#' # Compute PCA estimates
#' basis <- svd(example_data / sqrt(M), nu = M, nv = n)
#' initial_loadings <- sqrt(n) * basis$v[, 1:r]
#'
#' # Find minimum rotation using orthonormal basis initial_loadings
#' rotation_result <- find_local_factors(X = example_data, r = r, initial_loadings = initial_loadings)
#'
find_local_factors <- function(X, r, initial_loadings, parallel = FALSE, n_cores = NULL) {
  # Function to find the sparsest rotation of the leading Principal Components
  # of a M by n matrix X.
  # Under sparsity in the loading matrix this will identify the true loading matrix.
  #
  # returns three arguments:
  #   initial_loadings: Principal Component estimate
  #   rotated_loadings: Rotation of loading matrix with smallest l1-norm
  #   diagnostics: A list of diagnostics

  stopifnot(is.matrix(X) | is.data.frame(X))
  if("data.frame" %in% class(X)) X <- as.matrix(X)
  stopifnot(is.numeric(r))
  stopifnot(r %% 1 == 0 & r > 0)
  stopifnot(ncol(X) >= r)
  if(any(is.na(X)) | any(is.infinite(X))) stop("X cannot contain missing or infinite values.")
  if(!all(is.numeric(X))) stop("X cannot contain non-numeric values.")


  if(!missing(initial_loadings)){
    stopifnot(is.matrix(initial_loadings))
    stopifnot(ncol(initial_loadings) > 1)
    if(any(is.na(initial_loadings)) | any(is.infinite(initial_loadings))) stop("initial_loadings cannot contain missing or infinite values.")
    if(!all(is.numeric(initial_loadings))) stop("initial_loadings cannot contain non-numeric values.")
  }

  stopifnot(is.numeric(n_cores) | is.null(n_cores))
  if(is.numeric(n_cores)) stopifnot(n_cores %% 1 == 0)
  stopifnot(is.logical(parallel))
  if(parallel & is.null(n_cores)) stop("parallel set to TRUE but n_cores is NULL. Please specify n_cores for parallel execution.")
  if(!parallel & !is.null(n_cores)) warning("parallel set to FALSE but n_cores is not null. Defaulting to sequential execution.")

  n <- ncol(X)

  if (missing(initial_loadings)) {
    M <- nrow(X)
    svd_X <- svd(X / sqrt(M), nu = M, nv = n)
    eig_X <- svd_X$d^2
    initial_loadings <- sqrt(n) * svd_X$v[, 1:r]
  }

  if (any(round(t(initial_loadings) %*% initial_loadings) != diag(nrow = r)*n)){
    stop('The initial estimate initial_loadings should be an orthonormal basis of the loading space.
        Either drop argument (PCs will be used), or orthonormalize.')
  }

  # compute the rotated solution with minimal l1-norm
  rmat_min_results <- find_min_rotation(initial_loadings, parallel = parallel, n_cores = n_cores) #Finds solutions across a large grid of starting points
  rmat_min <- rmat_min_results$R
  rotation_results <- collate_solutions(rmat_min, initial_loadings, X) #Combine large number of solutions into candidates
  rotated_loadings <- rotation_results$rotated_loadings
  diagnostics <- rotation_results$diagnostics
  return(list(initial_loadings = initial_loadings, rotated_loadings = rotated_loadings, diagnostics = diagnostics))
}


