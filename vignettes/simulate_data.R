simulate_data_2d <- function(T, m1, m2, n, rho_t, rho_i, lambda_option, error_magnifier = 0.3) {
  # rho_t, rho_i correlation in error terms.
  
  # Dep defines the structure of the factors.
  Dep <- matrix(c(1, m1, n + 1 - m2, n), ncol = 2, byrow = TRUE)
  dim <- dim(Dep)
  r <- dim[1]
  F <- matrix(rnorm(T * r), nrow = T) %*% chol(matrix(c(1, 0.3, 0.3, 1), nrow = 2)) # Cholesky decomposition
  Lambda <- matrix(0, n, r)
  
  if (lambda_option == 1) {
    for (k in 1:r) {
      Lambda[Dep[k, 1]:Dep[k, 2], k] <- rep(1, Dep[k, 2] - Dep[k, 1] + 1) + rnorm(Dep[k, 2] - Dep[k, 1] + 1)
    }
  } else if (lambda_option == 2) {
    for (k in 1:r) {
      Lambda[Dep[k, 1]:Dep[k, 2], k] <- runif(Dep[k, 2] - Dep[k, 1] + 1, 0.1, 1.9)
    }
  } else if (lambda_option == 3) {
    for (k in 1:r) {
      Lambda[Dep[k, 1]:Dep[k, 2], k] <- rnorm(Dep[k, 2] - Dep[k, 1] + 1)
    }
  }
  
  # Create noise
  epsilon <- matrix(rnorm(T * n), nrow = T)
  cross_dep <- numeric(n)
  time_dep <- numeric(T)
  
  for (i in 1:n) {
    cross_dep[i] <- rho_i^(i - 1)
  }
  
  for (t in 1:T) {
    time_dep[t] <- rho_t^(t - 1)
  }
  
  A <- toeplitz(time_dep)
  B <- toeplitz(cross_dep)
  e <- sqrt(A) %*% epsilon %*% sqrt(B)
  
  # Now create X
  X <- F %*% t(Lambda) + error_magnifier * e
  
  # Also save normalized loading matrix
  truth_normal <- Lambda
  for (k in 1:r) {
    truth_normal[, k] <- sqrt(n) * Lambda[, k] / sqrt(sum(Lambda[, k]^2)) # normalize by columns to unit length
  }
  
  return(list(X = X, Lambda = Lambda, truth_normal = truth_normal, F = F, r = r))
}

# Constants -------
set.seed(200)

scaling=1
m1 = 1; 
m2 = 2; 
lambda_option=2; 
T <- 224
n <- 3
rho_t <- 0.3
rho_i <- 0.1
factorno <- 2

simulated_data <- simulate_data_2d(T, m1, m2, n, rho_t, rho_i, lambda_option, error_magnifier = 0.3)

lf <- l1rotation::local_factors(simulated_data$X, r = 2)

write.csv(simulated_data$X, "data.csv", row.names = FALSE)
write.csv(simulated_data$truth_normal, "truth_normal.csv", row.names = FALSE)
write.csv(lf$initial_loadings, "lambda0.csv", row.names = FALSE)



