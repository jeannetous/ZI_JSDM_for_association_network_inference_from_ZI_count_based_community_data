source("ZIPLN_generate_simulation_parameters.R")
library(MASS)

#' @description generates a named list with all the required parameters to simulate
#' data under the ZIPLN model
#' @param n number of rows in the Abundance matrix
#' @param p number of columns in the Abundance matrix
#' @param d number of covariates in the covariates matrix
#' @param omega_structure network structure for the precision matrix (erdos_renyi,
#' community or preferential_attachment)
#' @param sigma_sparse boolean, indicates whether the simulation should be based
#' on a sparse variance matrix - default is FALSE, corresponding to a sparse
#' precision matrix
#' @param zi_type type of zero-inflation
#' @param zi_covar_cluster boolean, if zi_type = covar, whether there should be
#' a division of the ZI values into rows and column clusters
#' @param min_X minimum value for X, either one single value for X, or a list of length d for each dimension
#' @param max_X maximum value for X, either one single value for X, or a list of length d for each dimension
#' @param mean_B mean of the Gaussian distributions B values are drawn from, can be a list with one mean per B dimension
#' @param sd_B standard deviation of the Gaussian distributions B values are drawn from, can be a list with one sd per B dimension
#' @param XB_max max value for exp(XB) mean of the Poisson distribution that Y follows
#' @param v calibration parameter to get Omega from a graph
#' @param u calibration parameter to get Omega from a graph
#' @param zi_mode_values if zi_type = sites or species, list of zi probabilities
#' @param proba_mode_zi list of probabilities of having each ZI contained in zi_mode_values
#' @param block_values if zi_type = "covar" and  zi_covar_cluster = TRUE, values
#' of X0 %*% B0 expected for each pair (row_cluster, col_cluster)
#' #' @param row_clusters if zi_type = "covar", divisions of 1:n into row
#' clusters, given as a list of labels, drawn at random if not specified
#' @param col_clusters if zi_type = "covar" divisions of 1:p into column
#' clusters, given as a list of labels, drawn at random if not specified
#' @param row_clusters_proba list of probabilities for row clusters, used only
#' if row_clusters=NULL, default is equiprobable distributions
#' @param col_clusters_proba list of probabilities for column clusters, used
#' only if col_clusters=NULL, default is equiprobable distributions
#' @param X0 optional, if zi_type = covar, list of ZI covariates
#' @param B0 optional, regression matrix for the ZI covariates, required if X0 is not null
#' @param min_X0 minimum value for X0, either one single value for X0, or a list of length d for each dimension,
#' applied only if zi_covar_cluster = FALSE (otherwise X0 is discrete)
#' @param max_X0 maximum value for X0, either one single value for X0, or a list of length d for each dimension
#' applied only if zi_covar_cluster = FALSE (otherwise X0 is discrete)
#' @param max_X0B0 maximum value for the mean of each column of X0 %*% B0
#' applied only if zi_covar_cluster = FALSE (otherwise X0 is discrete)
generate_all_ZIPLN_parameters <- function(n, p, d, add_intercept = TRUE,
                                          omega_structure = "erdos_renyi",
                                          sigma_sparse = FALSE,
                                          zi_type = c("covar", "sites", "species"),
                                          zi_covar_cluster = FALSE,
                                          min_X = 0, max_X = 10, mean_B  = 2,
                                          sd_B = 1, XB_max = 70,
                                          v = 0.3, u = 0.1,
                                          zi_mode_values = NULL,
                                          proba_mode_zi = NULL,
                                          block_values = NULL,
                                          row_clusters = NULL,
                                          col_clusters = NULL,
                                          row_clusters_proba = NULL,
                                          col_clusters_proba = NULL,
                                          X0 = NULL, B0 = NULL,
                                          min_X0 = 0, max_X0 = 10,
                                          max_X0B0 = -0.2){
  if(sigma_sparse){
    Sigma <- generate_sigma_sparse(p, omega_structure, v)
    Omega <- chol2inv(chol(Sigma))
  }else{
    Omega <- generate_omega(p, omega_structure, v, u)
    Sigma <- chol2inv(chol(Omega))
  }
  X <- generate_X(n, d, min_X, max_X, add_intercept)
  B <- generate_B(p, X, Sigma, mean_B, sd_B, XB_max)
  
  if(zi_type == "covar"){
    if(is.null(X0)){
      if(zi_covar_cluster){
        zi_params <- generate_X0_B0_cluster(n, p, block_values,
                                            row_clusters, col_clusters,
                                            row_clusters_proba, col_clusters_proba)
        B0 <- zi_params$B0 ; X0 <- zi_params$X0_num
        zi_params <- list(X0 = zi_params$X0, B0 = B0, X0_num = zi_params$X0_num)
      }else{
        X0 <- generate_X(n, d, min_X0, max_X0)
        colnames(X0) <- unlist(lapply(1:ncol(X0), f <- function(x) paste0("VZI", as.character(x))))
        B0 <- generate_B0(X0, max_X0B0)
        zi_params <- list(X0 = X0, B0 = B0)
      }
    }else{zi_params <- list(X0 = X0, B0 = B0)}
  }else{zi_params <- NULL}
  zi_proba <- generate_zi_proba(n, p, zi_type, zi_mode_values, proba_mode_zi,
                                X0, B0)
  return(list(Omega = Omega, Sigma = Sigma, X = X, B = B,
              zi_params = zi_params, zi_proba = zi_proba))
}

#' @description simulates data under the ZIPLN model for a list of fixed parameters
#' @param params named list of parameters for the ZIPLN model
simulate_ZIPLN_data <- function(params){
  n  <- nrow(params$X)
  mu <- params$X %*% params$B
  Y  <- rMLN(n, mu, Sigma = params$Sigma, N = rep(3000, n),
             zi_proba = params$zi_proba)
  return(Y)
}

#' Simulation of community count data
#'
#' Count data drawn under a possibly zero inflated multinomial model with a predefined (Gaussian) covariance structure.
#'
#' The workflow is:
#'  - draw latent abundances-basis under a centered multivariate normal with user-defined covariance matrix
#'  - use logistic transformation of abundance basis to proportions
#'  - draw counts (with a given sum N) under a multinomial distributions with previously defined proportions.
#'
#' @param n the sample size (number of communities to simulate)
#' @param mu vector of means of the latent variable (may containt covariates effect and offset); roughly equal to mean log-abundances
#' @param Sigma covariance matrix of the latent variables/species
#' @param N a vector of sequencing depth in each sample (default fixed to 3000 in all samples)
#' @param pi matrix of probabilities controling zero inflation in each sample (default fixed to 0 in all samples)
#' @return
#' @example
#' ## Simulation settings
#' p <- 50  # number fo communities
#' n <- 100 # sample size
#' d <- 3   # number of covariates
#' k <- p   # number of edges in the network
#' network <- sample_er(p, k) # true latent network (random structure)
#' covar_effect  <- 1
#' offset_effect <- 10  ## le larger, the less effect
#' Sigma <- graph2cov(network, v = 0.3, u = 0.1)
#' ## Covariate: One-way ANOVA with 3 modalities
#' group <- factor(sample(1:d, n, replace = TRUE))
#' X <- model.matrix(~ group + 0)
#' mu <- X %*% matrix(runif(d*p, -1,  1) , d, p)
#' ## sequencinf depth
#' N <- rnegbin(n, mu = 1000, theta = 2); N[N == 0] <- 1000
#' Y <- rMLN(n, mu = mu, Sigma = Sigma, N = N)
rMLN <- function(n, mu = matrix(0, n, p), Sigma, N = rep(3000, n),
                 zi_proba = rep(0, n)) {
  p <- ncol(Sigma)
  abundancies <- mu + mvrnorm(n, rep(0, p), as.matrix(Sigma))
  logistic <- function(x) { z <- exp(x); return(z / sum(z)) }
  proportions <- t(apply(abundancies, 1, logistic))
  ## sample counts
  ## first, convert matrix to list, for use with mapply
  prop.list <- lapply(seq_len(nrow(proportions)), function(i) proportions[i,])
  counts <- mapply(rmultinom, size = N, prob = prop.list, n = 1)
  W <- t(matrix(rbinom(n * p, size = 1, prob = zi_proba), n, p))
  counts <- counts * (W == 0)
  counts <- t(counts)
  
  which_zero  <- which(colSums(counts) == 0)
  if (length(which_zero) > 0) {
    counts[sample.int(n, 1), which_zero] <- 1
  }
  counts
}
