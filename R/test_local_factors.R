#' Test for the presence of local factors, as in [local_factors()], with additional flexibility.
#'
#' @inheritParams local_factors
#' @param loadings (optional) Matrix that represents a sparse basis of the loading space.
#'
#' @returns Returns a list with the following components:
#'  * `has_local_factors` Logical equal to `TRUE` if local factors are present.
#'  * `n_small` Integer denoting the number of small loadings in sparse rotation.
#'  * `gamma_n` Integer denoting the critical value to compare `n_small` to.
#'  * `h_n` Number denoting the cutoff used to determine which loadings are small.
#'  * `loadings` Matrix that is the rotation of the loadings that produces the smallest l1-norm (if not supplied).
#' @export
#'
#' @examples
#' # Minimal example with 2 factors, where X is a 224 by 207 matrix
#' r <- 2
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
#'    loadings = rotation_result$rotated_loadings
#' )
#'
#' test_result$has_local_factors
#'
test_local_factors <- function(X, r, loadings = NULL) {

  stopifnot(is.matrix(X) | is.data.frame(X))
  if("data.frame" %in% class(X)) X <- as.matrix(X)
  stopifnot(is.numeric(r))
  stopifnot(r %% 1 == 0 & r > 0)
  stopifnot(ncol(X) >= r)

  if(any(is.na(X)) | any(is.infinite(X))) stop("X cannot contain missing or infinite values.")
  if(!all(is.numeric(X))) stop("X cannot contain non-numeric values.")

  if(!is.null(loadings)){
    stopifnot(is.matrix(loadings))
    stopifnot(ncol(loadings) == r)
    stopifnot(ncol(loadings) > 1)
    if(any(is.na(loadings)) | any(is.infinite(loadings))) stop("Loadings matrix cannot contain missing or infinite values.")
    if(!all(is.numeric(loadings))) stop("Loadings matrix cannot contain non-numeric values.")
  } else  {
    loadings <- find_local_factors(X, r)$rotated_loadings
  }

  # Set hyperparameters
  n <- nrow(loadings)
  alpha_gamma <- 0.05
  c_gamma <- -1 * stats::qnorm(1 - alpha_gamma / 2, lower = FALSE)
  gamma0 <- 0.03
  h_n <- 1 / log(n)
  expected_small <- 1 / 2 * pracma::erfc(-h_n / sqrt(2)) - 1 / 2 * pracma::erfc(h_n / sqrt(2))
  gamma <- gamma0 + expected_small + c_gamma * sqrt((expected_small * (1 - expected_small)) / n)
  gamma_n <- floor(gamma * n)

  n_small <- colSums(abs(loadings) < h_n)
  most_small <- sort(n_small, decreasing = TRUE)
  has_local_factors <- most_small[1] > gamma_n


  return(list(has_local_factors = has_local_factors, n_small = n_small, gamma_n = gamma_n, h_n = h_n, loadings = loadings))
}


