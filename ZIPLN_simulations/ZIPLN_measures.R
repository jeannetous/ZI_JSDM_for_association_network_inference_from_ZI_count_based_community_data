library(Metrics)

#' @description computes partial correlations from precision
pcor_from_omega <- function(omega) {
  r <- -omega / tcrossprod(sqrt(diag(omega)))
  diag(r) <- 1
  r
}

#' @description computes recall, fallout, precision and f1-score for
#' an inferred network (precision matrices here) given the true one
#' @param omega_true true precision matrix
#' @param omega_estimate estimated precision matrix
roc_metrics <- function(omega_true, omega_estimate) {
  
  diag(omega_true) <- 0 ; p <- nrow(omega_true)
  roc <- function(theta) {
    diag(theta) <- 0
    
    nzero <- which(theta != 0)
    zero  <- which(theta == 0)
    
    true.nzero <- which(omega_true != 0)
    true.zero  <- which(omega_true == 0)
    
    TP <- sum(nzero %in% true.nzero)
    TN <- sum(zero %in%  true.zero) - p
    FP <- sum(nzero %in% true.zero)
    FN <- sum(zero %in%  true.nzero)
    recall    <- TP/(TP + FN) ## also recall and sensitivity
    fallout   <- FP/(FP + TN) ## also 1 - specificit
    precision <- TP/(TP + FP) ## also PPR
    f1_score  <- 2 * (precision * recall) / (precision + recall)
    recall[TP + FN == 0] <- NA
    fallout[TN + FP == 0] <- NA
    precision[TP + FP == 0] <- NA
    
    res <-  round(c(fallout,recall,precision, f1_score),3)
    res[is.nan(res)] <- 0
    names(res) <- c("fallout","recall", "precision", "f1_score")
    res
  }
  
  if (is.list(omega_estimate)) {
    return(as.data.frame(do.call(rbind, lapply(omega_estimate, roc))))
  } else {
    return(roc(omega_estimate))
  }
}

#' @description computes AUC from the roc measures
#' @param roc measures as computed by function roc_metrics, named list that contains
#' a list of fallout and a list of recall values
#' @param threshold from which the list of recall / fallout values should be cut
#' before only adding a (1, 1) point to the ROC curve
perf_auc <- function(roc, threshold = 1) {
  cut <- (roc$fallout < threshold) & (roc$recall < threshold)
  fallout <- c(0, roc$fallout[cut], threshold)
  recall  <- c(0, roc$recall[cut] , threshold)
  dx <- diff(fallout)
  res <- sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2
  res <- ifelse(is.character(res), NA, res)
  res
}

#' @description sparse precision path as obtained by the graphical-lasso, from
#' any variance-covariance matrix (including a singular one)
#' @param Sigma estimated covariance (can be low-rank)
#' @param n_penalties number of points along the "way" (number of graphical-lasso
#' penalties )
#' @param bound_by_spectrum boolean, if TRUE the sparse precision path is cut at 
#' Sigma's lowest non-zero eigen values:
#'  - Use FALSE for AUC computing: what matters is the ranking of the edges, not
#' the absolute values and setting it to TRUE can lead to a biased AUC. 
#'  - Use TRUE for point-wise estimates: when the penalty is too low, Omega
#'  estimation may "explode"
sparse_precision_path <- function(Sigma, n_penalties = 50,
                                  bound_by_spectrum = FALSE) {
  x  <- abs(Sigma[upper.tri(Sigma)])
  hi <- max(x)
  if (!is.finite(hi) || hi <= 0)
    return(rep(list(diag(1 / pmax(diag(Sigma), 1e-8))), n_penalties))
  lo <- max(min(x[x > 0]), hi * 1e-4)
  if (bound_by_spectrum) {
    ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    ev <- ev[ev > max(ev) * 1e-8]
    lo <- max(lo, min(ev))
  }
  lo <- min(lo, hi)
  rhos <- 10^seq(log10(hi), log10(lo), length.out = n_penalties)
  lapply(rhos, function(rho) glassoFast::glassoFast(Sigma, rho = rho)$wi)
}

#' @description selects a penalty along the sparsity way with the BIC, to give
#' a punctual estimate comparable with the BIC/StARS selection used for (ZI)PLN
#' @param Sigma estimated covariance
#' @param n number of observations (sites)
select_sparse_precision <- function(Sigma, n, n_penalties = 50) {
  path <- sparse_precision_path(Sigma, n_penalties, bound_by_spectrum = TRUE)
  bic  <- sapply(path, function(O) {
    if (anyNA(O)) return(NA_real_)
    ld <- determinant(O, logarithm = TRUE)
    if (ld$sign <= 0) return(NA_real_)          # Omega non definie positive
    k  <- sum(O[upper.tri(O)] != 0)
    -n * (as.numeric(ld$modulus) - sum(diag(Sigma %*% O))) + log(n) * k
  })
  if (all(!is.finite(bic)))                      # aucun point exploitable
    return(diag(1 / pmax(diag(Sigma), 1e-8)))    # repli : reseau vide
  path[[which.min(replace(bic, !is.finite(bic), Inf))]]
}

#' @description splits the error on a matrix (Sigma / Omega) into a
#' DIAGONAL part (variance errors) and an OFF-DIAGONAL part (network error, the
#' one this work focuses on)
matrix_rmse_parts <- function(matrix_hat, matrix_true) {
  d <- diag(TRUE, nrow(matrix_true))
  e <- (matrix_hat - matrix_true)^2
  c(matrix_rmse         = sqrt(mean(e)),
    matrix_rmse_diag    = sqrt(mean(e[d])),
    matrix_rmse_offdiag = sqrt(mean(e[!d])))
}

#' @description reference bound errors: error obtained by an estimator that 
#' infers an empty network (diagonal matrix).
empty_network_rmse <- function(matrix_true) {
  matrix_rmse_parts(diag(diag(matrix_true)), matrix_true)
}

#' @description given a precision matrix and a PLN model (collection of PLN fit
#' for different penalties), computes the AUC associated to the model
#' @param omega_true true precision matrix
#' @param model fitted PLN (with multiple penalties) or hmsc or gllvm model
#' @param omega_hat useful if hmsc or gllvm, the AUC is computed by thresholding
#' the values it contains
#' @param model_type whether model is a PLN, hmsc or gllvm model 
get_auc <- function(omega_true, model, omega_hat = NULL){
  roc <- roc_metrics(omega_true,
                     lapply(model$models, function(model) model$model_par$Omega))
  return(perf_auc(roc))
}

#' @description given a precision matrix and a low-rank variance-covariance matrix, 
#' computes the AUC associated to the low-rank variance-covariance matrix
#' @param omega_true true precision matrix
#' @param Sigma variance-covariance matrix (low-rank if gllvm / Hmsc)
get_auc_alternative <- function(omega_true, Sigma, n_penalties = 50){
  roc <- roc_metrics(omega_true, sparse_precision_path(Sigma, n_penalties))
  return(perf_auc(roc))
}

#' @description given a variance-covariance matrix and a model, 
#' computes the AUC 
#' @param sigma_true true variance-covariance matrix
#' @param model fitted hmsc model, if applicable
#' @param sigma_hat useful if PLN or gllvm, the AUC is computed by thresholding
#' the values it contains
#' @param model_type whether model is a PLN, hmsc or gllvm model 
get_sigma_auc <- function(sigma_true, model = NULL, sigma_hat = NULL,
                          model_type = c("PLN", "hmsc", "gllvm")){
  if(model_type == "hmsc"){
    pre_sigma  <- computeAssociations(model)
    roc <- roc_metrics(sigma_true,
                       lapply(seq(1, 0, length.out = 100),
                              function(posterior_proba){
                                sigma_hat  <- ((pre_sigma[[1]]$support > posterior_proba) + (pre_sigma[[1]]$support < (1-posterior_proba)) > 0) * pre_sigma[[1]]$mean
                                return(sigma_hat)}))
  }
  if(model_type %in% c("PLN", "gllvm")){
    sigma_quantiles <- quantile(sigma_hat[upper.tri(sigma_hat, diag = F)], probs = seq(1, 0, -0.05))
    roc <- roc_metrics(sigma_true,
                       lapply(sigma_quantiles,
                              function(q){
                                sigma_current <- sigma_hat
                                sigma_current[abs(sigma_current) <= q] <- 0
                                return(sigma_current)}))
  }
  return(perf_auc(roc))
}

#' @description plots the ROC curve plot (True Positive Rate as a function
#' of the False Positive Rate) of a model with different penalties
#' @param omega_true true precision matrix
#' @param PLN_model fitted PLN model (with multiple penalties)
plot_roc_curve <- function(omega_true, PLN_model){
  roc <- roc_metrics(omega_true,
                     lapply(PLN_model$models, function(model) model$model_par$Omega))
  plot(roc$recall, roc$fallout)
}

#' @description computes a collection of measures associated to a given model
#' @param model fitted PLN (with multiple penalties) / hmsc / gllvm model
#' @param params true parameters under which the data was simulated
#' @param unpenalised_model fitted ZIPLN model with no penalty
#' (useful if model_type = "ZIPLN")
#' @param model_selection for PLN - model selection criterion to compute the
#' measures for (BIC, ICL, StARS)
#' @param stability for PLN - when model_selection = StARS, level of stability
#' to use for stability selection
#' @param add_intercept whether there is an intercept, useful to get B rmse
#' @param AUC AUC if already known, to avoid recomputing it if it's already
#' been done with another model_selection criterion
#' @param AUC.alt AUC.alt if already known, to avoid recomputing it if it's
#' already been done with another model_selection criterion
#' @param Sigma_AUC Sigma_AUC if already known, to avoid recomputing it if it's
#' already been done with another model_selection criterion
#' @param model_type whether model is a PLN, hmsc or gllvm model 
get_measures <- function(model, params, unpenalised_model = NULL,
                         model_selection = NULL, stability = 0.8,
                         add_intercept = TRUE, AUC = NULL, AUC.alt = NULL,
                         Sigma_AUC = NULL, 
                         model_type = c("PLN", "ZIPLN", "hmsc", "gllvm")) {
  model_type <- match.arg(model_type)
  cor_true   <- cov2cor(params$Sigma)
  pcor_true  <- pcor_from_omega(params$Omega)
  
  # RMSE Upper bounds (common to all the models)
  omega_rmse_empty         = round(empty_network_rmse(params$Omega)[["matrix_rmse"]], 4)
  omega_rmse_empty_offdiag = round(empty_network_rmse(params$Omega)[["matrix_rmse_offdiag"]], 4)
  sigma_rmse_empty         = round(empty_network_rmse(params$Sigma)[["matrix_rmse"]], 4)
  sigma_rmse_empty_offdiag = round(empty_network_rmse(params$Sigma)[["matrix_rmse_offdiag"]], 4)
  pcor_rmse_empty          = round(empty_network_rmse(pcor_true)[["matrix_rmse"]], 4)
  pcor_rmse_empty_offdiag  = round(empty_network_rmse(pcor_true)[["matrix_rmse_offdiag"]], 4)
  cor_rmse_empty           = round(empty_network_rmse(cor_true)[["matrix_rmse"]], 4)
  cor_rmse_empty_offdiag   = round(empty_network_rmse(cor_true)[["matrix_rmse_offdiag"]], 4)
  
  if(model_type == "PLN" | model_type == "ZIPLN"){
    # Select best sparsity level according to the chosen criterion
    if(is.numeric(model_selection)){
      m <- model$getModel(model_selection)
    }else{
      if(model_selection == "StARS"){
        m <- model$getBestModel(model_selection, stability)
      }else{m <- model$getBestModel(model_selection)}
    }
    # B error
    if(add_intercept){
      PLN_B     <- m$model_par$B[c("(Intercept)", "V1"),]
    }else{PLN_B     <- m$model_par$B[c("V1"),]}
    
    if(nrow(PLN_B) == nrow(params$B)){
      B_rmse    <- Metrics::rmse(PLN_B, params$B)
    }else{B_rmse    <- Metrics::rmse(PLN_B, params$B[1:2,])}
    
    # Omega measures
    omega_hat          <- m$model_par$Omega
    omega_rmse_debiased <- NA_real_
    if(model_type == "PLN" | model_type == "ZIPLN"){
      if(model_type == "PLN"){
        S <- crossprod(m$var_par$M)/m$n + diag(colMeans(m$var_par$S**2), m$p, m$p)
      }else{
        S <- crossprod(m$var_par$M - params$X %*% m$model_par$B)/m$n + diag(colMeans(m$var_par$S * m$var_par$S))
      }
      omega_rmse_debiased <-
        matrix_rmse_parts(debias_network(S, m$model_par$Omega),
                          params$Omega)[["matrix_rmse"]]
      if(is.null(AUC)){
        AUC = get_auc(params$Omega, model)
        }
      if(is.null(AUC.alt) ){
        AUC.alt = get_auc_alternative(params$Omega,
                                      unpenalised_model$model_par$Sigma)
      }
      
      if(is.null(Sigma_AUC)){
        Sigma_AUC <- get_sigma_auc(params$Sigma,
                                   sigma_hat = unpenalised_model$model_par$Sigma ,
                                   model_type = "PLN")
      }
    }else{
      AUC     = get_auc_alternative(params$Omega, m$model_par$Sigma)
      AUC.alt = AUC
      Sigma_AUC <- get_sigma_auc(params$Sigma, model,
                                 sigma_hat = m$model_par$Sigma,
                                 model_type = "PLN")
    }
    
    ############################################################################
    
    omega_rmse_parts <- matrix_rmse_parts(omega_hat, params$Omega)
    omega_rmse       <- omega_rmse_parts[["matrix_rmse"]]
    pcor_rmse_parts  <- matrix_rmse_parts(pcor_from_omega(omega_hat), pcor_true)
    pcor_rmse        <- pcor_rmse_parts[["matrix_rmse"]]
    sigma_rmse_parts <- matrix_rmse_parts(m$model_par$Sigma, params$Sigma)
    sigma_rmse       <- sigma_rmse_parts[["matrix_rmse"]]
    cor_rmse_parts   <- matrix_rmse_parts(cov2cor(m$model_par$Sigma), cor_true)
    cor_rmse         <- cor_rmse_parts[["matrix_rmse"]]
    
    # Aggregated results
    res <- c(
      criterion = ifelse(model_selection == "StARS",
                         paste0(model_selection, "_", stability), model_selection),
      B_rmse = round(B_rmse, 4),
      omega_rmse = round(omega_rmse, 4),
      omega_rmse_diag    = round(omega_rmse_parts[["matrix_rmse_diag"]], 4),
      omega_rmse_offdiag = round(omega_rmse_parts[["matrix_rmse_offdiag"]], 4),
      omega_rmse_empty   = omega_rmse_empty,
      omega_rmse_empty_offdiag = omega_rmse_empty_offdiag,
      omega_rmse_debiased = round(omega_rmse_debiased, 4),
      pcor_rmse = round(pcor_rmse, 4),
      pcor_rmse_offdiag = round(pcor_rmse_parts[["matrix_rmse_offdiag"]], 4),
      pcor_rmse_empty   = pcor_rmse_empty,
      pcor_rmse_empty_offdiag = pcor_rmse_empty_offdiag,
      sigma_rmse = round(sigma_rmse, 4),
      sigma_rmse_offdiag = round(sigma_rmse_parts[["matrix_rmse_offdiag"]], 4),
      sigma_rmse_empty   = sigma_rmse_empty,
      sigma_rmse_empty_offdiag = sigma_rmse_empty_offdiag,
      cor_rmse = round(cor_rmse, 4),
      cor_rmse_offdiag = round(cor_rmse_parts[["matrix_rmse_offdiag"]], 4),
      cor_rmse_empty   = cor_rmse_empty,
      cor_rmse_empty_offdiag = cor_rmse_empty_offdiag,
      AUC = AUC, AUC.alt = AUC.alt,
      Sigma_AUC = Sigma_AUC,
      fit_rmse = rmse(m$fitted, params$Y),
      roc_metrics(params$Omega, omega_hat)
    )
  }
  
  if(model_type == "hmsc"){
    # B error
    if(add_intercept){
      hmsc_B     <- getPostEstimate(model, parName = "Beta")$mean[c(1,2), ]
    }else{hmsc_B     <- getPostEstimate(model, parName = "Beta")$mean[c(1), ]}
    
    if(nrow(hmsc_B) == nrow(params$B)){
      B_rmse    <- Metrics::rmse(hmsc_B, params$B)
    }else{B_rmse    <- Metrics::rmse(hmsc_B, params$B[1:2,])}
    
    # Omega measures
    pre_sigma           <- computeAssociations(model)
    sigma_glasso        <- pre_sigma[[1]]$mean
    omega_hat           <- select_sparse_precision(sigma_glasso, nrow(params$Y))
    omega_rmse_parts    <- matrix_rmse_parts(omega_hat, params$Omega)
    omega_rmse          <- omega_rmse_parts[["matrix_rmse"]]
    omega_rmse_debiased <- NA_real_
    pcor_rmse_parts     <- matrix_rmse_parts(pcor_from_omega(omega_hat), pcor_true)
    pcor_rmse           <- pcor_rmse_parts[["matrix_rmse"]]
    sigma_rmse_parts    <- matrix_rmse_parts(pre_sigma[[1]]$mean, params$Sigma)
    sigma_rmse          <- sigma_rmse_parts[["matrix_rmse"]]
    cor_rmse_parts      <- matrix_rmse_parts(cov2cor(pre_sigma[[1]]$mean), cor_true)
    cor_rmse            <- cor_rmse_parts[["matrix_rmse"]]
    
    if(is.null(AUC)){
      AUC     = get_auc_alternative(params$Omega, sigma_glasso)
      AUC.alt = NA_real_
    }
    Sigma_AUC <- get_sigma_auc(params$Sigma, model, 
                               model_type = "hmsc")
    
    # Fit measures
    fitted <- predict(model, XData = params$XData)
    fitted <- lapply(fitted,
                     function(Y){
                       Ybis <- Y
                       Ybis[is.na(params$Y_hurdle)] <- 0
                       return(Ybis)})
    fit_rmse <- min(unlist(lapply(fitted, function(Y) return(rmse(Y, params$Y)))))
    
    # Aggregated results
    res <- c(
      criterion = NA,
      B_rmse = round(B_rmse, 4),
      omega_rmse = round(omega_rmse, 4),
      omega_rmse_diag    = round(omega_rmse_parts[["matrix_rmse_diag"]], 4),
      omega_rmse_offdiag = round(omega_rmse_parts[["matrix_rmse_offdiag"]], 4),
      omega_rmse_empty   = omega_rmse_empty,
      omega_rmse_empty_offdiag = omega_rmse_empty_offdiag,
      omega_rmse_debiased = round(omega_rmse_debiased, 4),
      pcor_rmse = round(pcor_rmse, 4),
      pcor_rmse_offdiag = round(pcor_rmse_parts[["matrix_rmse_offdiag"]], 4),
      pcor_rmse_empty   = pcor_rmse_empty,
      pcor_rmse_empty_offdiag = pcor_rmse_empty_offdiag,
      sigma_rmse = round(sigma_rmse, 4),
      sigma_rmse_offdiag = round(sigma_rmse_parts[["matrix_rmse_offdiag"]], 4),
      sigma_rmse_empty   = sigma_rmse_empty,
      sigma_rmse_empty_offdiag = sigma_rmse_empty_offdiag,
      cor_rmse = round(cor_rmse, 4),
      cor_rmse_offdiag = round(cor_rmse_parts[["matrix_rmse_offdiag"]], 4),
      cor_rmse_empty   = cor_rmse_empty,
      cor_rmse_empty_offdiag = cor_rmse_empty_offdiag,
      AUC = AUC, AUC.alt = AUC.alt,
      Sigma_AUC = round(Sigma_AUC, 4),
      fit_rmse = round(fit_rmse, 4),
      roc_metrics(params$Omega, omega_hat)
    )
  }
  if(model_type == "gllvm"){
    
    # B error
    gllvm_B   <- model$params$Xcoef
    if(add_intercept){
      gllvm_B   <- model$params$Xcoef[, c("Intercept", "V1")]
    }else{gllvm_B   <- model$params$Xcoef[, c("V1")]}
    
    if(ncol(gllvm_B) == nrow(params$B)){
      B_rmse    <- Metrics::rmse(t(gllvm_B), params$B)
    }else{B_rmse    <- Metrics::rmse(t(gllvm_B), params$B[1:2,])}
    
    
    # Omega / Sigma measures
    Sigma               <- getResidualCov(model)$cov
    omega_hat           <- select_sparse_precision(Sigma, nrow(params$Y))
    omega_rmse_parts    <- matrix_rmse_parts(omega_hat, params$Omega)
    omega_rmse          <- omega_rmse_parts[["matrix_rmse"]]
    omega_rmse_debiased <- NA_real_
    pcor_rmse_parts     <- matrix_rmse_parts(pcor_from_omega(omega_hat), pcor_true)
    pcor_rmse           <- pcor_rmse_parts[["matrix_rmse"]]
    sigma_rmse_parts    <- matrix_rmse_parts(Sigma, params$Sigma)
    sigma_rmse          <- sigma_rmse_parts[["matrix_rmse"]]
    cor_rmse_parts      <- matrix_rmse_parts(cov2cor(Sigma), cor_true)
    cor_rmse            <- cor_rmse_parts[["matrix_rmse"]]
    
    if(is.null(AUC)){
      AUC     = get_auc_alternative(params$Omega, Sigma)
      AUC.alt = NA_real_
    }
    Sigma_AUC <- get_sigma_auc(params$Sigma, model, sigma_hat = Sigma,
                               model_type = "gllvm")
    
    # Fit measures
    fit_rmse   <- rmse(predict.gllvm(model), params$Y)
    
    # Aggregated results
    res <- c(
      criterion = NA,
      B_rmse = round(B_rmse, 4),
      omega_rmse = round(omega_rmse, 4),
      omega_rmse_diag    = round(omega_rmse_parts[["matrix_rmse_diag"]], 4),
      omega_rmse_offdiag = round(omega_rmse_parts[["matrix_rmse_offdiag"]], 4),
      omega_rmse_empty   = omega_rmse_empty,
      omega_rmse_empty_offdiag = omega_rmse_empty_offdiag,
      omega_rmse_debiased = round(omega_rmse_debiased, 4),
      pcor_rmse = round(pcor_rmse, 4),
      pcor_rmse_offdiag = round(pcor_rmse_parts[["matrix_rmse_offdiag"]], 4),
      pcor_rmse_empty   = pcor_rmse_empty,
      pcor_rmse_empty_offdiag = pcor_rmse_empty_offdiag,
      sigma_rmse = round(sigma_rmse, 4),
      sigma_rmse_offdiag = round(sigma_rmse_parts[["matrix_rmse_offdiag"]], 4),
      sigma_rmse_empty   = sigma_rmse_empty,
      sigma_rmse_empty_offdiag = sigma_rmse_empty_offdiag,
      cor_rmse = round(cor_rmse, 4),
      cor_rmse_offdiag = round(cor_rmse_parts[["matrix_rmse_offdiag"]], 4),
      cor_rmse_empty   = cor_rmse_empty,
      cor_rmse_empty_offdiag = cor_rmse_empty_offdiag,
      AUC = AUC, AUC.alt = AUC.alt,
      Sigma_AUC = Sigma_AUC,
      fit_rmse = round(fit_rmse, 4),
      roc_metrics(params$Omega, omega_hat)
    )
  }
  res
}