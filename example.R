pacman::p_load(
  readr,
  dplyr,
  devtools
)

load_all()


X <- read_csv("../example_data/example_data1.csv", col_names = FALSE) %>%
  as.matrix()

ex1 <- local_factors(X, 4)

X <- read_csv("../example_data/example_data2.csv", col_names = FALSE) %>%
  as.matrix()

ex2 <- local_factors(X, 8)


