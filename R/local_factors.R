#' Test whether local factors are present in a given dataset `X` and return the rotation of the loading matrix with the smallest l1-norm.
#'
#' @param X A (usually standardized) t by n matrix of observations.
#' @param r An integer denoting the number of factors in X.
#'
#' @returns Returns a list with the following components:
#'  * `has_local_factors` a logical equal to `TRUE` if local factors are present
#'  * `Lambda0` the principal component estimate of the loading matrix
#'  * `Lambda_rotated` a matrix that is the rotation of the loading matrix that produces the smallest l1-norm.
#'  * `rotation_diagnostics` a list containing 3 components"
#'      * `R` the rotation matrix that when used to rotate `Lambda0` produces the smallest l1-norm.
#'      * `l1_norm` a vector of length `r` containing the value of the l1 norm each solution generates
#'      * `sol_frequency` a vector of length `r` containing the frequency in the initial grid of each solution
#'
#' @export
#'
local_factors <- function(X, r) {

  #set.seed(1909)
  T <- nrow(X)
  n <- ncol(X)

  # Compute PCA estimates
  pca <- svd(X / sqrt(T))
  eig_X <- pca$d^2
  Lambda0 <- sqrt(n) * pca$v[, 1:r]

  # Test for local factors and estimate rotated loading matrix
  result <- test_local_factors(X, r)
  has_local_factors <- result$has_local_factors
  Lambda_rotated <- result$Lambda

  # Illustrate loading matrices
  pc_plot <- plot_loading_matrix(Lambda0, xlab = "k", title = "Principal Component estimate")
  pc_rotated_plot <- plot_loading_matrix(Lambda_rotated, xlab = "k", title = "Rotated estimate (l1-criterion)")
  small_loadings_plot <- plot_small_loadings(result)

  return(list(
    has_local_factors = has_local_factors,
    Lambda0 = Lambda0,
    Lambda_rotated = Lambda_rotated,
    rotation_diagnostics = result$rotation_diagnostics,
    pc_plot = pc_plot,
    pc_rotated_plot = pc_rotated_plot,
    small_loadings_plot = small_loadings_plot))
}

plot_loading_matrix <- function(data, xlab = "", ylab = "", title = ""){

  if(is.matrix(data)) data <- convert_mat_to_df(data) %>% glimpse()

  scale_fill_pal <- select_palette(data$value, type = "level")

  ggplot2::ggplot(data, ggplot2::aes(column, as.numeric(row), fill = value)) +
    ggplot2::geom_tile() +
    scale_fill_pal +
    ggplot2::labs(x = xlab, y = ylab, title = title, fill = "") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      barwidth = 1,
      barheight = 15,
      label.theme = ggplot2::element_text(size = 12),
      draw.ulim = TRUE,
      draw.llim = TRUE)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

plot_small_loadings <- function(result, xlab = "k", ylab = "", title = ""){

  n_small <- result$n_small
  gamma <- result$gamma
  h_n <- result$h_n
  n_cols <- ncol(result$Lambda)

  tibble::as_tibble(n_small) %>%
    dplyr::rownames_to_column(var = "factor") %>%
    dplyr::mutate(factor = as.numeric(factor)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = factor, y = value), size = 3) +
    ggplot2::geom_hline(yintercept = gamma, linetype = "dashed", size = 1) +
    ggplot2::ylim(c(min(gamma - 10, min(n_small - 5)), max(gamma + 5, max(n_small) + 5))) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::xlim(c(1, n_cols)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}


convert_mat_to_df <- function(mat){

  df <- mat %>% tibble::as_tibble() %>% dplyr::rownames_to_column(var = "row") %>%
    dplyr::mutate(row = factor(row, levels = 1:nrow(mat))) %>%
    tidyr::pivot_longer(starts_with("V"), names_to = "column") %>%
    dplyr::mutate(column = factor(stringr::str_remove(column, "V"), levels = 1:ncol(mat)))

  return(df)
}



select_palette <- function(range, type, breaks = NULL){

  # Define midpoint
  mp <- switch(type,
               difference = 0,
               ratio = 1,
               level = median(range)
  )

  min <- min(range, na.rm = TRUE)
  max <- max(range, na.rm = TRUE)

  # Get limits
  d <- max(abs(mp - min), abs(max - mp))
  if (min < mp & max > mp) {
    auto_limits <- c(mp - d, mp + d)
    colors <- c("#008EFF", "#EBF0F4", "#961046")
  }
  if (min >= mp){
    auto_limits <- c(mp, max)
    colors <- c("#EBF0F4", "#961046")

  }
  if (max < mp) {
    auto_limits <- c(mp - d, mp)
    colors <- c("#008EFF", "#EBF0F4")

  }

  palette <- ggplot2::scale_fill_gradient2(
    high = scales::muted("darkblue"), mid = "white", low = "maroon",
    limits = auto_limits, midpoint = mp, oob = scales::squish
  )
  return(palette)
}



