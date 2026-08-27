library(igraph)

####################### Functions to generate Omega / Sigma ####################

#' @description generates an Erdos-Renyi graph
#' @param p number of nodes in the graph
#' @param prob probability of having an edge between 2 nodes
erdos_renyi_graph <- function(p, prob = 0.05){
  as_adjacency_matrix(sample_gnp(p, prob))
}

#' @description generates a preferential attachment graph
preferential_attachment_graph <- function(p){
  as_adjacency_matrix(sample_pa(p, m = 1, directed = FALSE))
}

#' @description generates a community structure attachment graph
#' @param prob probability of belonging to each group in the community
#' @param prob_in probability of having an edge between 2 nodes from the same group
#' @param prob_in probability of having an edge between 2 nodes from different groups
community_graph <- function(p, prob = c(1/2,1/4,1/4), prob_in = 0.1, prob_out = 0.05) {
  pref_mat <- matrix(prob_out, length(prob), length(prob))
  diag(pref_mat) <- prob_in
  graph_mat <- as_adjacency_matrix(sample_sbm(p,
                                              pref.matrix = pref_mat,
                                              block.sizes = c(rmultinom(1, p, prob)) ))
  graph_mat
}

#' @description generates a sparse precision matrix Omega
#' @param p dimension of the graph
#' @param omega_structure type of structure for the underlying graph
#' @param v calibration parameter to get Omega from a graph
#' @param u calibration parameter to get Omega from a graph
generate_omega <- function(p, omega_structure, v = 0.3, u = 0.1){
  cond <- FALSE
  while(!cond){
    if(omega_structure == "erdos_renyi") G <- erdos_renyi_graph(p)
    if(omega_structure == "preferential_attachment") G <- preferential_attachment_graph(p)
    if(omega_structure == "community") G <- community_graph(p)
    
    # Ensuring that the network is not empty for AUC to make sense
    if(max(G) == 0){
      off_diag_indices <- which(row(matrix(1:p, p, p)) != col(matrix(1:p, p, p)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      G[selected_index[["row"]], selected_index[["col"]]] <- 1
      G[selected_index[["col"]], selected_index[["row"]]] <- 1
    }
    omega_tilde <- G * v
    omega <- omega_tilde + diag(abs(min(eigen(omega_tilde)$values)) + u, p, p)
    
    # Ensuring that the network is not full for AUC to make sense
    if(min(omega) > 0){
      off_diag_indices <- which(row(matrix(1:p, p, p)) != col(matrix(1:p, p, p)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      omega[selected_index[["row"]], selected_index[["col"]]] <- 0
      omega[selected_index[["col"]], selected_index[["row"]]] <- 0
    }
    
    # Including some variability + negative values in omega
    upper <- upper.tri(omega, diag = FALSE)
    pos_upper <- omega[upper] > 0
    upper_pos <- matrix(FALSE, nrow = p, ncol = p)
    upper_pos[upper] <- pos_upper
    to_replace <- rnorm(sum(upper_pos), mean = omega[upper_pos], sd = 0.2)
    omega[upper_pos] <- unlist(lapply(1:length(to_replace), f <- function(i){
      ifelse(to_replace[[i]] > 0.9 | to_replace[[i]] < 0,
             omega[upper_pos][[i]],
             to_replace[[i]])}))
    
    omega[upper_pos] <- round(omega[upper_pos], 2)
    omega <- as.matrix(omega)
    prob <- runif(sum(upper_pos))
    to_negate <- prob < 0.4
    omega[upper_pos][to_negate] <- - omega[upper_pos][to_negate]
    omega[lower.tri(omega)] <- t(omega)[lower.tri(omega)]
    
    # controlling for omega to be positive definite and its inverse to not have too high values
    cond <- ! is.complex(eigen(omega)$values )
    if(cond) cond <- all(eigen(omega)$values > 0)
    if(cond) cond <- max(abs(solve(omega))) < 3
  }
  as.matrix(omega)
}

#' @description generates a sparse variance-covariance matrix Sigma
#' @param p dimension of the graph
#' @param sigma_structure type of structure for the underlying graph
#' @param v calibration parameter to get Sigma from a graph
generate_sigma_sparse <- function(p, sigma_structure, v = 0.3){
  repeat {
    if(sigma_structure == "erdos_renyi") G <- erdos_renyi_graph(p)
    if(sigma_structure == "preferential_attachment") G <- preferential_attachment_graph(p)
    if(sigma_structure == "community") G <- community_graph(p)
    
    # Ensuring that the network is not empty for AUC to make sense
    if(max(G) == 0){
      off_diag_indices <- which(row(matrix(1:p, p, p)) != col(matrix(1:p, p, p)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      G[selected_index[["row"]], selected_index[["col"]]] <- 1
      G[selected_index[["col"]], selected_index[["row"]]] <- 1
    }
    G <- as.matrix(G)
    S <- G * v
    S[upper.tri(S)] <- S[upper.tri(S)] * ifelse(runif(sum(upper.tri(S))) < .4, -1, 1)
    S[lower.tri(S)] <- t(S)[lower.tri(S)]
    # S <- S + diag(abs(min(eigen(S)$values)) + 0.4, p)
    S <- S + diag(abs(min(eigen(S)$values)) + 0.4 + rnorm(p, 0.5, 0.1), p)
    if (all(eigen(S)$values > 0)) return(S)
  }
}
################## Functions to add zero-inflation in the data #################
#' @description generate B0 that will define ZI probabilities with X0
#' @param X0 matrix of ZI-related covariates
#' @param max_X0B0 maximum value for the mean of each column of X0 %*% B0
generate_B0 <- function(X0, max_X0B0 = -0.2){
  d <- ncol(X0)
  B0 <- matrix(rep(1, d*p), nrow=d)
  for(dim in 1:d){B0[dim,] = runif(p, min=-1, max = 1)}
  correcting_factors <- unlist(lapply(colMeans(X0 %*% B0),
                                      f <- function(x){ifelse(x <= max_X0B0, 1,
                                                              max_X0B0 / x)}))
  B0 <- sweep(B0, 2, correcting_factors, `*`)
}

#' @description generates a discrete X0 and the corresponding B0 to have zi-proba
#' defined by rows and columns clusters
#' @param n number of rows in the final zi_proba matrix
#' @param block_values values of X0 %*% B0 expected for each pair (row_cluster, col_cluster)
#' @param row_clusters divisions of 1:n into row clusters, given as a list of
#' labels, drawn at random if not specified
#' @param col_clusters divisions of 1:p into column clusters, given as a list of
#' labels, drawn at random if not specified
#' @param row_clusters_proba list of probabilities for row clusters, used only if row_clusters=NULL,
#' default is equiprobable distributions
#' @param col_clusters_proba list of probabilities for column clusters, used only if col_clusters=NULL,
#' default is equiprobable distributions
generate_X0_B0_cluster <- function(n, p, block_values,
                                   row_clusters = NULL, col_clusters = NULL,
                                   row_clusters_proba = NULL,
                                   col_clusters_proba = NULL) {
  a <- nrow(block_values) ; b <- ncol(block_values)
  if(is.null(row_clusters)){
    if(is.null(row_clusters_proba)){
      row_clusters <- sort(rep(1:a, length.out = n))
    }else{row_clusters <- sort(sample(1:a, size = n, replace = TRUE, prob = row_clusters_proba))}
  }
  if(is.null(col_clusters)){
    if(is.null(col_clusters_proba)){
      col_clusters <- sort(rep(1:b, length.out = p))
    }else{col_clusters <- sort(sample(1:b, size = p, replace = TRUE, prob = col_clusters_proba))}
  }
  B0 <- block_values[, col_clusters, drop = FALSE]
  B0 <- apply(B0, c(1,2), f <- function(x){rnorm(1, x, 0.05)})
  X0 <- generate_discrete_X(n, 1, a, row_clusters)
  colnames(X0) <- unlist(lapply(1:ncol(X0), f <- function(x) paste0("VZI", as.character(x))))
  X0_num <- model.matrix(~ . - 1, data = as.data.frame(X0))
  return(list("X0" = X0, "X0_num" = X0_num, "B0" = B0))
}

#' @description generates a matrix of zero-inflation probabilities
#' @param n number of rows in the output matrix
#' @param p number of columns in the output matrix
#' @param zi_type zero-inflation type (covariate-dependent, site-wise = row-wise or species-wise = column-wise)
#' @param zi_mode_values if zi_type = sites or species, list of zi probabilities
#' @param proba_mode_zi list of probabilities of having each ZI contained in zi_mode_values
#' @param X0 if zi_type = covar, list of ZI covariates
#' @param B0 if zi_type = covar, regression parameters, so that zi_proba = logit(X0 %*% B0)
generate_zi_proba <- function(n, p, zi_type = c("covar", "sites", "species"),
                              zi_mode_values = NULL, proba_mode_zi = NULL,
                              X0 = NULL, B0 = NULL){
  zi_type <- match.arg(zi_type)
  if(zi_type == "covar"){
    X0B0 <- X0 %*% B0
    zi_proba <- exp(X0B0) / (1 + exp(X0B0))
  }else{
    n_mode_zi_proba <- length(zi_mode_values)
    if(is.null(proba_mode_zi)){
      if(zi_type == "sites"){groups <- sort(rep(1:n_mode_zi_proba, length.out = n))
      }else{groups <- sort(rep(1:n_mode_zi_proba, length.out = p))}
    }else{
      if(zi_type == "sites"){groups <- sample(1:n_mode_zi_proba, size = n, replace = TRUE, prob = proba_mode_zi)
      }else{groups <- sample(1:n_mode_zi_proba, size = p, replace = TRUE, prob = proba_mode_zi)}
    }
    if(zi_type == "sites"){
      zi_proba_list <-  unlist(lapply(1:n, f <- function(i){rnorm(1, mean = zi_mode_values[[groups[[i]]]], sd = 0.05)}))
      zi_proba <- matrix(rep(zi_proba_list, p), nrow = n, byrow = F)
    }else{
      zi_proba_list <-  unlist(lapply(1:p, f <- function(j){rnorm(1, mean = zi_mode_values[[groups[[j]]]], sd = 0.05)}))
      zi_proba <- matrix(rep(zi_proba_list, n), nrow = n, byrow = T)
    }
  }
  zi_proba <- apply(zi_proba, c(1, 2), f <- function(x) min(0.95, max(0.05, x)))
  return(zi_proba)
}

#' @description adds 0s to an abundance matrix
#' @param Y abundance matrix
#' @param matrix of zero-inflation probabilities for each (row, column) in Y
add_zero_inflation <- function(Y, zi_proba){
  Z <- apply(zi_proba, c(1, 2), f <- function(x) rbinom(1,1,x))
  Y[Z == 1] <- 0
  return(Y)
}

##################### Functions to generate other parameters ###################

#' @description generates continuous covariates matrix X
#' @param n number of rows in X
#' @param d number of column in X
#' @param min_X minimum value for X, either one single value for X, or a list of length d for each dimension
#' @param max_X maximum value for X, either one single value for X, or a list of length d for each dimension
generate_X <- function(n, d, min_X = 0, max_X = 1, add_intercept = TRUE){
  if(length(min_X == 1)) min_X <- rep(min_X, d)
  if(length(max_X == 1)) max_X <- rep(max_X, d)
  X = matrix(rep(1, n * d), nrow=n)
  for(dim in 1:d){X[,dim] = runif(n, min=min_X[[dim]], max = max_X[[dim]])}
  if(add_intercept){
    X <- cbind(rep(1, n), X)
    colnames(X) <- c("Intercept", unlist(lapply(1:d, f <- function(x) paste0("V", as.character(x)))))
  }else{colnames(X) <- unlist(lapply(1:d, f <- function(x) paste0("V", as.character(x))))}
  return(X)
}

#' @description generates discrete covariates matrix X
#' @param n number of rows in X
#' @param d number of column in X
#' @param n_cat_values number of values to include, either one single value for X,
#' or a list of length d for each dimension, the values are distributed equiprobably in X
#' @param row_clusters clustering of the rows (NULL by default - then the clusters
#' are drawn at random)
generate_discrete_X <- function(n, d, n_cat_values, row_clusters = NULL){
  X = matrix(rep(1, n * d), nrow=n)
  if(is.null(row_clusters)){
    for(dim in 1:d){X[,dim] = sort(rep(1:n_cat_values[[dim]], length.out = n))}
  }else{for(dim in 1:d){X[,dim] = row_clusters}}
  X <- apply(X, c(1, 2), f <- function(x) LETTERS[as.numeric(x)])
  return(X)
}

#' @description generates regression matrix for a given covariate matrix
#' @param p number of columns in B
#' @param X covariates matric
#' @param Sigma variance-covariance matrix used in the model
#' @param SNR signal to noise ratio, ratio between Sigma's variance and that of XB
generate_B <- function(p, X, Sigma, mean_B  = 2, sd_B = 1, XB_max = 70){
  d <- ncol(X)
  B <- matrix(rep(1, d*p), nrow=d)
  for(dim in 1:d){B[dim,] = rnorm(p, mean = mean_B, sd = sd_B)}
  B <- (log(XB_max) / max(X %*% B)) * B
  return(B)
}
