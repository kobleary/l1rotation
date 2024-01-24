
normalize <- function(X, p = 2){
  stopifnot(is.matrix(X))
  norms <- apply(X, p = p, 2, pracma::Norm)
  X_norm <- sweep(X, 2, norms, FUN = "/")
  return(X_norm)
}

# Returns the norm of each column of a matrix
vecnorm <- function(X, p = 2){
  stopifnot(is.matrix(X))
  apply(X, p = p, 2, pracma::Norm)
}


find_min_rotation <- function(Lambda) {
  r <- ncol(Lambda)
  no_draws <- gridsize(r)
  fval <- rep(0, no_draws)
  exitflag <- rep(0, no_draws)

  # Create starting points for algorithm
  initial_draws <- matrix(rnorm(r * no_draws), nrow = r)
  initial_draws <- normalize(initial_draws, p = 2)

  # Convert to polar coordinates (OA.7)?
  theta <- matrix(0, nrow = r - 1, ncol = no_draws)
  #for (kk in 1:(r - 2)) {
  #  theta[kk, ] <- acot(initial_draws[kk, ] / vecnorm(initial_draws[(kk + 1):r, ]))
  #}
  #theta[r - 1, ] <- 2 * acot((initial_draws[(r - 1), ] + vecnorm(initial_draws[(r - 1):r, ])) / initial_draws[r, ])

  # RK: Using Wikipedia definition

  theta <- cartesian_to_spherical(initial_draws)

  #theta <- matrix(0, nrow = r - 1, ncol = no_draws)
  #for (kk in 1:(r - 2)) {
  #  theta[kk, ] <- atan2( vecnorm(initial_draws[(kk + 1):r, ]), initial_draws[kk, ])
  #}
  #theta[r - 1, ] <- atan2( initial_draws[r, ], initial_draws[(r - 1), ] )


  # Optimization in polar coordinates happens w.r.t. theta
  angles <- theta
  l <- nrow(angles)
  for (rep in 1:no_draws) {
    starting_point <- theta[, rep]
    result <- optim(
      starting_point,
      objectivefcn_spherical, Lambda = Lambda,
      control = list(maxit = 200 * l),
      method = 'Nelder-Mead'
      )
    angles[, rep] <- result$par
    fval[rep] <- result$value
    exitflag[rep] <- result$convergence
  }

  # Convert back to cartesian coordinates, need to edit to generalize across minimum number of factors (requires at least 2)
  R <- spherical_to_cartesian(angles)

  #R <- matrix(0, nrow = r, ncol = no_draws)
  #R[1, ] <- cos(angles[1, ])
  #for (kk in 2:(r - 1)) {
  #  R[kk, ] <- col_prod(sin(angles[1:(kk - 1), ]))*cos(angles[kk, ])
  #}
  #R[r, ] <- col_prod(sin(angles))

  return(list(R = R, fval = fval, exitflag = exitflag))
}

col_prod <- function(data){
  if(is.matrix(data)) matrixStats::colProds(data)
  else{
    c(data)
  }
}

cartesian_to_spherical <- function(X){
  r <- nrow(X)
  no_draws <- ncol(X)

  theta <- matrix(0, nrow = r - 1, ncol = no_draws)
  for (kk in 1:(r - 2)) {
    theta[kk, ] <- atan2( vecnorm(X[(kk + 1):r, ]), X[kk, ])
  }
  theta[r - 1, ] <- atan2( X[r, ], X[(r - 1), ] )

  return(theta)
}

spherical_to_cartesian <- function(theta){
  r <- nrow(theta) + 1
  no_draws <- ncol(theta)

  R <- matrix(0, nrow = r, ncol = no_draws)

  R[1, ] <- cos(theta[1, ])

  if(r == 2){
    R[2, ] <- sin(theta[1, ]) * cos(theta[2, ])
  }
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
  n <- nrow(Lambda)
  r <- ncol(Lambda)
  R <- rep(0, r)
  R[1] <- cos(theta[1])
  for (kk in 2:(r-1)) {
    R[kk] <- prod(sin(theta[1:(kk-1)])) * cos(theta[kk])
  }
  R[r] <- prod(sin(theta))
  f <- sum(abs(Lambda %*% R))
  return(f)
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
