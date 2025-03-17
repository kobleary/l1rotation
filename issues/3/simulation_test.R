# Simulation test
library(tidyverse)
library(matrixStats)

devtools::load_all()

data_path <- here::here("issues", "3", "datastore", "simulation")
out_path <- here::here("issues","3")
r_tbl_approx_sparsity <- tibble()
matlab_tbl_approx_sparsity <- tibble()
sims <- 100:500
r_tbl <- tibble()
matlab_tbl <- tibble()

for(i in sims){

  X <- read_csv(file.path(data_path, str_glue("X_{i}.csv")), col_names = FALSE) %>%
    as.matrix()

  true_lambda <- read_csv(file.path(data_path, str_glue("true_lambda_{i}.csv")), col_names = FALSE) %>%
    as.matrix()

  lambda_matlab <- read_csv(file.path(data_path, str_glue("lambda_rotated_{i}.csv")), col_names = FALSE) %>%
    as.matrix()

  lambda_rotated <- find_local_factors(X, 4)$Lambda_rotated

  r_tbl <- bind_rows(
    r_tbl,
    colMaxs(abs(t(lambda_rotated) %*% true_lambda)/nrow(true_lambda))
  )

  matlab_tbl <- bind_rows(
    matlab_tbl,
    colMaxs(abs(t(lambda_matlab) %*% true_lambda)/nrow(true_lambda))
  )

}

for(i in sims){

  X <- read_csv(file.path(data_path, str_glue("X_approx_spars_{i}.csv")), col_names = FALSE) %>%
    as.matrix()

  true_lambda <- read_csv(file.path(data_path, str_glue("true_lambda_approx_spars_{i}.csv")), col_names = FALSE) %>%
    as.matrix()

  lambda_matlab <- read_csv(file.path(data_path, str_glue("lambda_rotated_approx_spars_{i}.csv")), col_names = FALSE) %>%
    as.matrix()

  lambda_rotated <- find_local_factors(X, 4)$Lambda_rotated

  r_tbl_approx_sparsity <- bind_rows(
    r_tbl_approx_sparsity,
    colMaxs(abs(t(lambda_rotated) %*% true_lambda)/nrow(true_lambda))
  )

  matlab_tbl_approx_sparsity <- bind_rows(
    matlab_tbl_approx_sparsity,
    colMaxs(abs(t(lambda_matlab) %*% true_lambda)/nrow(true_lambda))
  )

}

# compile results into data frame
tbl <- bind_rows(
  r_tbl %>% mutate(implementation = "R"),
  matlab_tbl %>% mutate(implementation = "Matlab")
) %>%
  pivot_longer(cols = -implementation, names_to = "factor", names_prefix = "X", values_to = "MC")

tbl_approx_spars <- bind_rows(
  r_tbl_approx_sparsity %>% mutate(implementation = "R"),
  matlab_tbl_approx_sparsity %>% mutate(implementation = "Matlab")
) %>%
  pivot_longer(cols = -implementation, names_to = "factor", names_prefix = "X", values_to = "MC")


write_csv(tbl, here::here(out_path, "r_matlab_sim_boxplot.csv"))
write_csv(tbl_approx_spars, here::here(out_path, "r_matlab_approx_spars_sim_boxplot.csv"))


# Plot and save boxplot
theme_set(theme_minimal())

boxplot <- tbl %>%
  ggplot(aes(factor, MC, fill = implementation, color = implementation)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch = FALSE, alpha = 0.5) +
  labs(color = "", fill = "", x = "Factor", y = expression(MC[l])) +
  theme(legend.position = c(.8, .2)) +
  scale_color_manual(values = c("R" = "cornflowerblue", "Matlab" = "darkorange")) +
  scale_fill_manual(values = c("R" =" cornflowerblue", "Matlab" = "darkorange"))

boxplot

ggsave(plot = boxplot, filename = here::here(out_path, "r_matlab_sim_boxplot.png"))


boxplot_approx <- tbl_approx_spars %>%
  ggplot(aes(factor, MC, fill = implementation, color = implementation)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch = FALSE, alpha = 0.5) +
  labs(color = "", fill = "", x = "Factor", y = expression(MC[l])) +
  theme(legend.position = c(.8, .2)) +
  scale_color_manual(values = c("R" = "cornflowerblue", "Matlab" = "darkorange")) +
  scale_fill_manual(values = c("R" =" cornflowerblue", "Matlab" = "darkorange"))

boxplot_approx

ggsave(plot = boxplot_approx, filename = here::here(out_path, "r_matlab_sim_approx_boxplot.png"))

