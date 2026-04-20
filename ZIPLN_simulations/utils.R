sigmoid <- function(x){ 1 / (1 + exp(-x))}

logit <- function(p){ log(p/(1 - p))}


debias_network <- function(S, network) {
  zeros <- which((network == 0 & upper.tri(network)), arr.ind = T)
  if(nrow(zeros) != 0) {
    network_debias <- glasso::glasso(S, rho = 0.0001, penalize.diagonal = FALSE, zero = zeros)$wi
  } else {
    network_debias <- network
  }
  return(network_debias)
}
