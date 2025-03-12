#' Test for the presence of local factors, as in [local_factors()], with additional flexibility.
#'
#' @inheritParams local_factors
#' @param Lambda (optional) Matrix that represents a sparse basis of the loading space.
#' @param alpha_gamma (optional) A numeric value that represents the tuning parameter (default = 0.05).
#'
#' @returns Returns a list with the following components:
#'  * `has_local_factors` Logical equal to `TRUE` if local factors are present.
#'  * `n_small` Integer denoting the number of small loadings in sparse rotation.
#'  * `gamma_n` Integer denoting the critical value to compare `n_small` to.
#'  * `h_n` Number denoting the cutoff used to determine which loadings are small.
#'  * `Lambda` Matrix that is the rotation of the loading matrix that produces the smallest l1-norm.
#' @export
#'
#' @examples
#' # Minimal example with 4 factors, where X is a 500 by 300 matrix
#' r <- 4
#' M <- nrow(example_data)
#' n <- ncol(example_data)
#'
#' # Find minimum rotation
#' rotation_result <- find_local_factors(X = example_data, r)
#'
#' # Test if sparse basis has local factors
#' test_result <- test_local_factors(
#'    X = example_data,
#'    r = r,
#'    Lambda = rotation_result$Lambda_rotated,
#'    alpha_gamma = 0.05
#' )
#'
#' test_result$has_local_factors
#'
test_local_factors <- function(X, r, Lambda = NULL, alpha_gamma = 0.05) {

  stopifnot(is.matrix(X) | missing(X))
  stopifnot(is.numeric(r))
  if(r %% 1 != 0) stop(stringr::str_glue("Specified r = {r} is not an integer."))
  stopifnot(is.numeric(alpha_gamma))

  if(!is.null(Lambda)){

    stopifnot(ncol(Lambda) == r)

    n <- nrow(Lambda)
    if (any(round(diag(t(Lambda) %*% Lambda)) != rep(n, r))) {
      stop('The initial estimate Lambda0 should be an orthonormal basis of the loading space.
        Either drop argument (PCs will be used), or orthonormalize.') }
  } else  {
    rotation_results <- find_local_factors(X, r)
    rotation_diagnostics <- rotation_results$diagnostics
    Lambda <- rotation_results$Lambda_rotated
    n <- nrow(Lambda)
  }

  c_gamma <- -1 * stats::qnorm(1 - alpha_gamma / 2, lower = FALSE)
  gamma0 <- 0.03
  h_n <- 1 / log(n)
  expected_small <- 1 / 2 * pracma::erfc(-h_n / sqrt(2)) - 1 / 2 * pracma::erfc(h_n / sqrt(2))
  gamma <- gamma0 + expected_small + c_gamma * sqrt((expected_small * (1 - expected_small)) / n)
  gamma_n <- floor(gamma * n)
  n_small <- colSums(abs(Lambda) < h_n)
  most_small <- sort(n_small, decreasing = TRUE)
  has_local_factors <- most_small[1] > gamma_n


  return(list(has_local_factors = has_local_factors, n_small = n_small, gamma_n = gamma_n, h_n = h_n, Lambda = Lambda))
}


