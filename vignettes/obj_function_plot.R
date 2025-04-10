r <- 2
no_draws <- 2000
Lambda0 <- read.csv('https://raw.githubusercontent.com/kobleary/quarto-example/refs/heads/main/data/lambda0.csv') %>%
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


# # slider_input <- seq(-3, 3, by = .1)
# #
# # obj_y <- map_dbl(slider_input, \(value) objectivefcn_spherical(value, lf$Lambda0))
#
# obj_plot <- obj_function_data %>%
#   ggplot() +
#   geom_point(aes(theta_orig, y_orig), color = "orange", size = 3) +
#   geom_line(aes(theta, y)) +
#   theme_minimal() +
#   labs(y = "", title = "Objective function")


