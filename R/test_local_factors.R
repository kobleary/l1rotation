#' Test for local factors, as in `local_factors` with additional flexibility.
#'
#' @inheritParams local_factors
#' @param Lambda (optional) a matrix that represents a sparse basis of the loading space
#' @param eig_X (optional) a vector of eigenvalues of X'X/t in decreasing order (only used if X is empty)
#' @param alpha_gamma (optional) a numeric value (default = 0.05)
#'
#' #' @returns Returns a list with the following components:
#'  * `has_local_factors` a logical equal to `TRUE` if local factors are present
#'  * `n_small` an integer denoting the number of small loadings in sparse rotation
#'  * `gamma_n` an integer denoting the critical value to compare `n_small` to.
#'  * `h_n` a number denoting the cutoffused to determine which loadings are small
#'  * `Lambda` rotation of PCs with smallest l1-norm
#'  * `rotation_diagnostics` a list containing 3 components"
#'      * `R` the rotation matrix that when used to rotate `Lambda0` produces the smallest l1-norm.
#'      * `l1_norm` a vector of length `r` containing the value of the l1 norm each solution generates
#'      * `sol_frequency` a vector of length `r` containing the frequency in the initial grid of each solution
#'
#' @export
test_local_factors <- function(X, r, Lambda = NULL, eig_X = NULL, alpha_gamma = 0.05) {
  # Function to test whether X has local factors
  #
  # returns up to four arguments:
  #    has_local_factors: Logical equal to one if local factors are present
  #    n_small: Number of small loadings in sparse rotation
  #    gamma_n: Critical value to compare n_small to.
  #    h_n: Threshold used for "small" loadings
  #    Lambda: Rotation of loading matrix with smallest l1-norm
  # See README.txt for more detail

  ## Preliminaries
  T <- nrow(X)
  n <- ncol(X)

  if (is.null(Lambda)) {
    rotation_results <- find_local_factors(X, r)
    rotation_diagnostics <- rotation_results$diagnostics
    Lambda <- rotation_results$Lambda_rotated

  }

  if (is.null(X)) {
    n <- nrow(Lambda)
    if (any(round(diag(t(Lambda) %*% Lambda)) != rep(r, n))) {
      stop("Loading matrix may not be properly normalized. Consider only passing two arguments (X,r).")
    }
    if (is.null(eig_X)) {
      stop("If no data X is supplied, at least eigenvalues needed to determine critical values")
    }
  } else {
    eig_X <- sort(eigen(t(X) %*% X / T)$values, decreasing = TRUE)
  }

  if (is.null(alpha_gamma)) {
    alpha_gamma <- 0.05
  }

  c_gamma <- -1 * qnorm(1 - alpha_gamma / 2, lower = FALSE)
  gamma0 <- 0.03
  h_n <- 1 / log(n)
  expected_small <- 1 / 2 * erfc(-h_n / sqrt(2)) - 1 / 2 * erfc(h_n / sqrt(2))
  gamma <- gamma0 + expected_small + c_gamma * sqrt((expected_small * (1 - expected_small)) / n)
  gamma_n <- floor(gamma * n)
  n_small <- colSums(abs(Lambda) < h_n)
  most_small <- sort(n_small, decreasing = TRUE)
  has_local_factors <- most_small[1] > gamma_n

  return(list(has_local_factors = has_local_factors, n_small = n_small, gamma_n = gamma_n, h_n = h_n, Lambda = Lambda, rotation_diagnostics = rotation_diagnostics))
}


