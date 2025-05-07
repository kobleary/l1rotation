library(dplyr)
r <- 2
no_draws <- 2000
Lambda0 <- read.csv('https://raw.githubusercontent.com/kobleary/l1rotation/refs/heads/main/vignettes/lambda0.csv') %>%
  as.matrix()

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


write.csv(obj_function_data, "obj_function_data.csv")

