load_matrix <- function(path){
  read.csv(path, header = FALSE) |>
    as.matrix()
}
