find_min_rotation <- function(Lambda) {
  stopifnot(is.matrix(Lambda))
  stopifnot(ncol(Lambda) > 1)

  r <- ncol(Lambda)
  no_draws <- gridsize(r)
  l1_norm <- rep(0, no_draws)
  exitflag <- rep(0, no_draws)

  # Create starting points for algorithm
  initial_draws <- matrix(stats::rnorm(r * no_draws), nrow = r)
  initial_draws <- normalize(initial_draws, p = 2)

  theta <- cartesian_to_spherical(initial_draws)

  # Optimization in polar coordinates happens w.r.t. theta
  angles <- theta
  l <- nrow(angles)
  for (rep in 1:no_draws) {
    starting_point <- theta[, rep]
    result <- stats::optim(
      starting_point,
      objectivefcn_spherical, Lambda = Lambda,
      control = list(maxit = 200 * l),
      method = 'Nelder-Mead'
      )
    angles[, rep] <- result$par
    l1_norm[rep] <- result$value
    print(result$value)
    exitflag[rep] <- result$convergence
  }

  # Convert back to cartesian coordinates, need to edit to generalize across minimum number of factors (requires at least 2)
  R <- spherical_to_cartesian(angles)

  return(list(R = R, l1_norm = l1_norm, exitflag = exitflag))
}


normalize <- function(X, p = 2){
  stopifnot(is.matrix(X))
  norms <- apply(X, p = p, 2, pracma::Norm)
  X_norm <- sweep(X, 2, norms, FUN = "/")
  return(X_norm)
}

# Returns the norm of each column of a matrix
vecnorm <- function(X, p = 2){
  #stopifnot(is.matrix(X))
  if(!is.matrix(X)){
    pracma::Norm(X[(kk):r, ], p = 2)
  }
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
  #stopifnot(all(theta >= 0))
  print(theta)
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

objectivefcn_spherical <- function(theta, Lambda) {
  R <- spherical_to_cartesian(theta)
  sum(abs(Lambda %*% R))
}

gridsize <- function(factorno) {
  # defines number of random draws to start search for local minma from
  if (factorno == 2) {
    no_randomgrid <- 300
  } else if (factorno == 3) {
    no_randomgrid <- 500
  } else if (factorno == 4) {
    no_randomgrid <- 1000
  } else if (factorno == 5) {
    no_randomgrid <- 2000
  } else if (factorno > 5 && factorno < 9) {
    no_randomgrid <- 3000
  } else {
    no_randomgrid <- 5000
  }

  return(no_randomgrid)
}
