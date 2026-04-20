source("ZIPLN_simulate_data.R")

simulate_negative_binomial_data <- function(params){
  n      <- nrow(params$X) ; p <- nrow(params$Sigma)
  Z      <- mvrnorm(n, rep(0, p), as.matrix(params$Sigma))
  mu     <- params$X %*% params$B + Z
  lambda <- exp(mu)
  U      <- matrix(rgamma(n * p, params$r, params$r), nrow = n)
  counts  <- apply(lambda * U, c(1,2), function(x) rpois(1, x))
  
  # Adding zero-inflation
  W <- matrix(rbinom(n * p, size = 1, prob = params$zi_proba), n, p)
  counts <- counts * (W == 0)
  
  which_zero  <- which(colSums(counts) == 0)
  if (length(which_zero) > 0) {
    counts[sample.int(n, 1), which_zero] <- 1
  }
  
  return(counts)
}