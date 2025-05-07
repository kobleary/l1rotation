# Note: must be run from "vignettes" folder

library(dplyr)
source("helper_functions.R")

# Simple simulation -------
set.seed(200)
T <- 224
n <- 3
r <- 2
rho <- .4
primitive <- matrix(rnorm(T * 2), nrow = T, ncol = 2)

factor_cov <- matrix(c(1, rho, 0, sqrt(1-rho^2)), byrow = TRUE, nrow = 2)
F <- primitive %*% factor_cov
Lambda <- matrix(c(1, .5, 1, 0, 0, -1), ncol = 2, byrow = TRUE)
e <- matrix(rnorm(T * n), nrow = T, ncol = n)
X <- F %*% t(Lambda) + 0.5 * e

lf <- l1rotation::local_factors(X, r = 2)


r <- 2
no_draws <- 2000
Lambda0 <- lf$initial_loadings

# Construct objective function ---------
# Create starting points for algorithm
initial_draws <- matrix(stats::rnorm(r * no_draws), nrow = r)
initial_draws <- normalize(initial_draws, p = 2)

theta <- cartesian_to_spherical(initial_draws)

slider_input <- seq(-3.13, 3.13, by = .01)


y <- purrr::map_dbl(
  slider_input,
  \(x) objectivefcn_spherical(x, Lambda0)
)


obj_function_data <- tibble::tibble(theta = slider_input, y = y) %>%
  mutate(
    dist1 = abs(theta - round(2*pi,2)),
    dist2 = abs(theta - round(pi, 2)/2),
  ) %>%
  mutate(min_dist1 = dist1 == min(dist1),
         min_dist2 = dist2 == min(dist2),
         highlight = min_dist1 | min_dist2
  )

obj_function_data %>%
  ggplot(aes(theta, y)) + geom_line() +
  geom_point(
    data = filter(obj_function_data, highlight == TRUE),
    aes(theta, y)
    )

# Write all to csv in "vignettes" folder ----
write.csv(X, "data.csv", row.names = FALSE)
write.csv(Lambda, "truth_normal.csv", row.names = FALSE)
write.csv(lf$initial_loadings, "lambda0.csv", row.names = FALSE)
write.csv(obj_function_data, "obj_function_data.csv")




