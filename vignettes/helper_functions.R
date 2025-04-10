# find_min_rotation functions ------

normalize <- function(X, p = 2){
  stopifnot(is.matrix(X))
  norms <- apply(X, p = p, 2, pracma::Norm)
  X_norm <- sweep(X, 2, norms, FUN = "/")
  return(X_norm)
}

# Returns the norm of each column of a matrix
vecnorm <- function(X, p = 2){
  apply(X, p = p, 2, pracma::Norm)
}


col_prod <- function(data){
  if(is.matrix(data)) matrixStats::colProds(data)
  else{
    c(data)
  }
}

# Assumes radius is equal to 1 (that is, X is normalized)
cartesian_to_spherical <- function(X){
  stopifnot(nrow(X) > 1)
  r <- nrow(X)
  no_draws <- ncol(X)

  theta <- matrix(0, nrow = r - 1, ncol = no_draws)
  if(r-2 > 0){
    for (kk in 1:(r - 2)) {
      theta[kk, ] <- atan2( vecnorm(X[(kk + 1):r, ]), X[kk, ])
    }
  }
  theta[r - 1, ] <- atan2( X[r, ], X[(r - 1), ] )

  return(theta)
}


spherical_to_cartesian <- function(theta){
  if(!is.matrix(theta)) {
    r <- length(theta) + 1
    R <- rep(0, r)
    R[1] <- cos(theta[1])
    if(r > 2){
      for (kk in 2:(r-1)) {
        R[kk] <- prod(sin(theta[1:(kk-1)])) * cos(theta[kk])
      }
    }
    R[r] <- prod(sin(theta))
    return(R)
  }
  stopifnot(nrow(theta) > 0)

  r <- nrow(theta) + 1
  no_draws <- ncol(theta)

  R <- matrix(0, nrow = r, ncol = no_draws)

  R[1, ] <- cos(theta[1, ])

  if(r > 2){
    for (kk in 2:(r - 1)) {
      R[kk, ] <- col_prod(sin(theta[1:(kk - 1), ]))*cos(theta[kk, ])

    }
  }

  if(r > 1){
    R[r, ] <- col_prod(sin(theta))
  }

  return(R)
}

objectivefcn_spherical <- function(theta, initial_loadings) {
  R <- spherical_to_cartesian(theta)
  sum(abs(initial_loadings %*% R))
}




# Plotly Plotting functions -----------------------------------------------
# Plotting functions ---------

add_points <- function(p, data, color = "black", size = 2){
  return(
    add_trace(p,
              type = "scatter3d",
              x = data[,1],
              y = data[,2],
              z = data[,3],
              mode = "markers",
              marker = list(
                size = size,
                color = color
              ),
              showlegend = FALSE
    )
  )
}

add_plane <- function(p, spanning_vectors, color = "blue", opacity = 0.5){

  origins <- matrix(c(
    0, 0, 0,
    0, 0, 0
  ), ncol = 3, byrow = TRUE)

  v1 <- spanning_vectors[1,]
  v2 <- spanning_vectors[2,]

  # Create parameter grid
  scale = 2
  t <- seq(-scale, scale, length.out = 20)
  s <- seq(-scale, scale, length.out = 20)

  # Calculate plane points
  origin <- c(0, 0, 0)
  plane_points <- expand.grid(t = t, s = s) %>%
    mutate(
      x = origin[1] + v1[1] * t + v2[1] * s,
      y = origin[2] + v1[2] * t + v2[2] * s,
      z = origin[3] + v1[3] * t + v2[3] * s
    )

  # Reshape for surface plot
  X <- matrix(plane_points$x, nrow = length(t))
  Y <- matrix(plane_points$y, nrow = length(t))
  Z <- matrix(plane_points$z, nrow = length(t))

  return(add_surface(p,
                     x = X,
                     y = Y,
                     z = Z,
                     colorscale = list(c(0, color), c(1, color)),
                     showscale = FALSE,
                     opacity = opacity
  ))

}


create_arrow_head <- function(origin, vector) {

  arrow_scale = 0.1
  arrow_angle = 30

  # Calculate vector endpoint
  endpoint <- origin + vector

  # Normalize vector
  vector_length <- sqrt(sum(vector^2))
  unit_vector <- vector / vector_length

  # Create perpendicular vectors for arrow head
  if (!all(unit_vector == c(0, 0, 1))) {
    perp1 <- pracma::cross(unit_vector, c(0, 0, 1))
  } else {
    perp1 <- pracma::cross(unit_vector, c(1, 0, 0))
  }
  perp1 <- perp1 / sqrt(sum(perp1^2))
  perp2 <- pracma::cross(unit_vector, perp1)

  # Arrow parameters
  arrow_length <- vector_length * arrow_scale
  arrow_angle_rad <- arrow_angle * pi / 180

  # Create arrow head points
  arrow_point1 <- endpoint - arrow_length * (
    unit_vector * cos(arrow_angle_rad) +
      perp1 * sin(arrow_angle_rad)
  )
  arrow_point2 <- endpoint - arrow_length * (
    unit_vector * cos(arrow_angle_rad) -
      perp1 * sin(arrow_angle_rad)
  )

  rbind(arrow_point1, endpoint, arrow_point2)
}



add_vectors <- function(p, vectors, origins = NULL, width = 3, color = "gray"){

  if(is.null(origins)){
    origins <- matrix(c(
      0, 0, 0,
      0, 0, 0
    ), ncol = 3, byrow = TRUE)
  }

  for (i in 1:nrow(origins)) {
    # Add vector line
    p <- p %>%
      add_trace(
        type = "scatter3d",
        x = c(origins[i,1], origins[i,1] + vectors[i,1]),
        y = c(origins[i,2], origins[i,2] + vectors[i,2]),
        z = c(origins[i,3], origins[i,3] + vectors[i,3]),
        mode = "lines",
        line = list(color = color, width = width),
        showlegend = FALSE
      )

    # Add arrow head
    arrow_points <- create_arrow_head(origins[i,], vectors[i,])
    p <- p %>%
      add_trace(
        p,
        type = "scatter3d",
        x = arrow_points[,1],
        y = arrow_points[,2],
        z = arrow_points[,3],
        mode = "lines",
        line = list(color = color, width = width),
        showlegend = FALSE
      )
  }

  return(p)
}

make_base_plot <- function(data, vectors, point_size = 2, vector_width = 7, vector_color = "gray"){

  plot_3d <- plot_ly() %>%
    layout(
      scene = list(
        camera = list(
          eye = list(x = .7, y = .7, z = .7)
        ),
        xaxis = list(title = "X", showgrid = FALSE),
        yaxis = list(title = "Y", showgrid = FALSE),
        zaxis = list(title = "Z", showgrid = FALSE),
        aspectmode = "cube"
      ),
      margin = list(l = 0, r = 0, t = 0, b = 0)
    ) %>%
    add_points(data = data, size = point_size) %>%
    add_plane(spanning_vectors = vectors) %>%
    add_vectors(vectors = vectors, width = vector_width, color = vector_color)

  return(plot_3d)

}
