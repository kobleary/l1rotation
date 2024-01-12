collate_solutions <- function(rmat_min, Lambda_0, eig_X) {
  n <- nrow(Lambda_0)
  factorno <- nrow(rmat_min)
  no_randomgrid <- ncol(rmat_min)
  epsilon_rot <- 0.05
  
  l1_min <- colSums(abs(Lambda_0 %*% rmat_min))
  sort_index <- order(l1_min)
  l1_min_sort <- l1_min[sort_index]
  rmat_min_sort <- rmat_min[, sort_index]
  rmat_min_sort <- rmat_min_sort * repmat(sign(rmat_min_sort[1, ]), factorno, 1) # same as in matlab up until here
  
  tic() # This step takes 44 seconds with 8 factors and 3000 draws (probably can be improved)
  norms <- matrix(0, nrow = no_randomgrid, ncol = no_randomgrid)
  for (i in 1:no_randomgrid) {
    for (j in i:no_randomgrid) {
      
      norms[i,j] <- (pracma::Norm(rmat_min_sort[, i] - rmat_min_sort[, j], p = 2) / factorno)
      if (norms[i,j] < epsilon_rot) {
      #if ((pracma::Norm(rmat_min_sort[, i] - rmat_min_sort[, j], p = 2) / factorno) < epsilon_rot) {
        rmat_min_sort[, j] <- rmat_min_sort[, i]
        l1_min_sort[j] <- l1_min_sort[i]
      }
    }
  }
  toc()
  
  candidate_tibble <- rmat_min_sort %>% 
    t() %>% as.data.frame() %>% 
    mutate(l1_norm = l1_min_sort) %>%
    group_by(across(everything())) %>%  
    dplyr::count() %>% 
    arrange(l1_norm) %>% 
    ungroup()
  
  
  #l1_min <- colSums(abs(Lambda_0 %*% rmat_min)) # calculate l1 norms that the rotations found in rmat_min generate
  
  # candidate_tibble <- rmat_min %>% t() %>% as.data.frame() %>% 
  #   mutate(l1_norm = l1_min) %>% 
  #   arrange(l1_norm)
  
  # rmat_min_sort <- candidate_tibble %>% select(-l1_norm) %>% as.matrix() %>% t() 
  # rmat_min_sort <- rmat_min_sort * repmat(sign(rmat_min_sort[1, ]), factorno, 1) # Normalize sign
  # 
  # candidate_tibble_ <- rmat_min_sort %>% t() %>% as.data.frame() %>% cbind(candidate_tibble$l1_norm)
  # 
  # for (i in 1:no_randomgrid) {
  #   for (j in i:no_randomgrid) {
  #     if (pracma::Norm(rmat_min_sort[, i] - rmat_min_sort[, j], p = 2) / factorno < epsilon_rot) {
  #       rmat_min_sort[, j] <- rmat_min_sort[, i]
  #       l1_min_sort[j] <- l1_min_sort[i]
  #     }
  #   }
  # }
  
  # starting_indices <- match(l1_min_sort, unique(l1_min_sort))
  # l1_min_unique <- l1_min_sort[starting_indices]
  # group_id <- rep(1:length(starting_indices), times = diff(c(starting_indices, length(l1_min_sort) + 1)))
  # value_counts <- cbind(unique(l1_min_unique), tabulate(group_id))
  # 
  # rmat_min_unique <- rmat_min_sort[, starting_indices]
  # 
  
  # Construct R from the candidates
  rmat_min_unique <- candidate_tibble %>% dplyr::select(-c(l1_norm, n)) %>% as.matrix() %>% t()
  Lambda_rotated <- Lambda_0 %*% rmat_min_unique[, 1]
  R <- rmat_min_unique[, 1]
  #sol_frequency <- value_counts[1, 2]
  #fval <- value_counts[1, 1]
  
  for (kk in 2:nrow(candidate_tibble)) {
    temp <- cbind(Lambda_rotated, Lambda_0 %*% rmat_min_unique[, kk])
    not_singular <- min(eigen(t(temp) %*% temp)$values) / n > sqrt(1 / factorno) / 3 && 
      min(eigen(t(temp) %*% temp)$values) / n > min(eigen(t(Lambda_rotated) %*% Lambda_rotated)$values) / n / 4
    if(not_singular){
      Lambda_rotated <- temp
      R <- cbind(R, rmat_min_unique[, kk])
     # sol_frequency <- c(sol_frequency, value_counts[kk, 2])
    #  fval <- c(fval, value_counts[kk, 1])
    }
  }

  candidate_tibble_prune <- R %>% t() %>% as.data.frame() %>%
    left_join(candidate_tibble)
  
  r_temp <- ncol(Lambda_rotated)
  I <- diag(factorno)
  l1_norm_update <- c()
  while (r_temp < factorno) {
    min_eig <- rep(0, factorno)
    for (ell in 1:factorno) {
      temp <- cbind(Lambda_rotated, Lambda_0[, ell])
      min_eig[ell] <- min(eigen(t(temp) %*% temp)$values)
    }
    index <- which.min(min_eig * sqrt(eig_X[1:factorno]))
    print(index)
    print(min_eig)
    Lambda_rotated <- cbind(Lambda_rotated, Lambda_0[, index])
    R <- cbind(R, I[, index])
    l1_norm_update <- c(l1_norm_update, sum(abs(Lambda_0[, index])))
    r_temp <- ncol(Lambda_rotated)
  }
  
  candidate_tibble_prune <- R %>% t() %>% as.data.frame() %>%
    left_join(candidate_tibble) %>%
    mutate(
      l1_norm = c(l1_norm[!is.na(l1_norm)], l1_norm_update),
      n = ifelse(is.na(n), 0, n)
    ) %>% glimpse()
  
  #return(list(R = R, fval = fval, sol_frequency = sol_frequency, Lambda_rotated = Lambda_rotated))
  return(list(R = R, fval = candidate_tibble_prune$l1_norm, sol_frequency = candidate_tibble_prune$n, Lambda_rotated = Lambda_rotated))
  
}


