library(Matrix)
library(colorspace)
library(ggplot2)
library(viridis)
library(scales)

walk(
  c(
    "test_local_factors.R",
    "find_local_factors.R",
    "find_min_rotation.R",
    "collate_solutions.R"
  ),
  source
)

local_factors <- function(X, r) {
  
  set.seed(1909)
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
  
  ggplot(data, aes(column, as.numeric(row), fill = value)) +
    geom_tile() +
    #scale_fill_viridis(name = "", option = "viridis") +
    scale_fill_pal +
    labs(x = xlab, y = ylab, title = title, fill = "") +
    guides(fill = guide_colorbar(
      barwidth = 1,
      barheight = 15,
      label.theme = element_text(size = 12),
      draw.ulim = TRUE,
      draw.llim = TRUE)) +
    theme_minimal() +
    theme(
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) 
}

plot_small_loadings <- function(result, xlab = "k", ylab = "", title = ""){
  
  n_small <- result$n_small
  gamma <- result$gamma
  h_n <- result$h_n
  n_cols <- ncol(result$Lambda)

  as_tibble(n_small) %>%
    rownames_to_column(var = "factor") %>% 
    mutate(factor = as.numeric(factor)) %>% 
    ggplot() +
    geom_point(aes(x = factor, y = value), size = 3) +
    geom_hline(yintercept = gamma, linetype = "dashed", size = 1) +
    ylim(c(min(gamma - 10, min(n_small - 5)), max(gamma + 5, max(n_small) + 5))) +
    labs(x = xlab, y = ylab, title = title) +
    xlim(c(1, n_cols)) +
    theme_minimal() +
    theme(legend.position = "none")
}


convert_mat_to_df <- function(mat){
  
  df <- mat %>% as_tibble() %>% rownames_to_column(var = "row") %>% 
    mutate(row = factor(row, levels = 1:nrow(mat))) %>% 
    pivot_longer(starts_with("V"), names_to = "column") %>% 
    mutate(column = factor(str_remove(column, "V"), levels = 1:ncol(mat)))
  
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
  
  palette <- scale_fill_gradient2(
    high = muted("darkblue"), mid = "white", low = "maroon", 
    limits = auto_limits, midpoint = mp, oob = squish
  )
  
  #eturn(list(midpoint = mp, limits = auto_limits, colors = colors))
  return(palette)
}



