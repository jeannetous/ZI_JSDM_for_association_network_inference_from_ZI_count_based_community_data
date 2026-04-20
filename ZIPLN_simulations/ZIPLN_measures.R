library(Metrics)

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


#' @description given a precision matrix and a PLN model (collection of PLN fit
#' for different penalties), computes the AUC associated to the model
#' @param omega_true true precision matrix
#' @param model fitted PLN model (with multiple penalties) or fitted hmsc model
#' @param omega_hat useful if hmsc or gllvm, the AUC is computed by thresholding
#' the values it contains
#' @param model_type whether model is a PLN or hmsc model 
get_auc <- function(omega_true, model, omega_hat = NULL,
                    model_type = c("PLN", "hmsc", "gllvm")){
  model_type <- match.arg(model_type)
  if(model_type == "PLN"){
    roc <- roc_metrics(omega_true,
                       lapply(model$models, function(model) model$model_par$Omega))
  }
  # if(model_type == "hmsc"){
  #   pre_Sigma = computeAssociations(model)
  #   roc <- roc_metrics(omega_true,
  #                      lapply(0.1*(1:10) - 0.05,
  #                             function(posterior_proba){
  #                               Sigma <- ((pre_Sigma[[1]]$support > posterior_proba) + (pre_Sigma[[1]]$support < (1-posterior_proba)) > 0) * pre_Sigma[[1]]$mean
  #                               Omega <- solve(Sigma)
  #                               return(Omega)}))
  # }
  if(model_type %in% c("hmsc", "gllvm")){
    omega_quantiles <- quantile(omega_hat[upper.tri(omega_hat, diag = F)], probs = seq(0, 1, 0.05))
    roc <- roc_metrics(omega_true,
                       lapply(omega_quantiles,
                              function(q){
                                omega_current <- omega_hat
                                omega_current[abs(omega_current) <= q] <- 0
                                return(omega_current)}))
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

#' @description computes a collection of measures associated to a given PLN model
#' @param model fitted PLN model (with multiple penalties) / fitted hmsc model
#' @param params true parameters under which the data was simulated
#' @param model_selection for PLN - model selection criterion to compute the measures for (BIC, ICL, StARS)
#' @param stability for PLN - when model_selection = StARS, level of stability to use for stability selection
#' @param posterior_proba for hmsc - posterior proba to select Sigma non-zero values
#' @param add_intercept whether there is an intercept, useful to get B rmse
#' @param AUC AUC if already known, to avoid recomputing it if it's already been done with another model_selection criterion
#' @param model_type whether model is a PLN, hmsc or gllvm model 
get_measures <- function(model, params, model_selection = NULL,
                         stability = 0.8, posterior_proba = 0.9,
                         add_intercept = TRUE, AUC = NULL,
                         model_type = c("PLN", "ZIPLN", "hmsc", "gllvm")) {
  model_type <- match.arg(model_type)
  
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
    ### Debiasing omega estimate ###############################################
    if(model_type == "PLN"){
      S <- crossprod(m$var_par$M)/m$n + diag(colMeans(m$var_par$S**2), m$p, m$p)
    }else{
      S <- crossprod(m$var_par$M - params$X %*% m$model_par$B)/m$n + diag(colMeans(m$var_par$S * m$var_par$S))
    }
    omega_hat <- debias_network(S, m$model_par$Omega)
    ############################################################################

    omega_rmse <- Metrics::rmse(omega_hat, params$Omega)
    if(is.null(AUC)) AUC = get_auc(params$Omega, model, model_type)
    
    # Aggregated results
    res <- c(
      criterion = ifelse(model_selection == "StARS",
                         paste0(model_selection, "_", stability), model_selection),
      B_rmse = round(B_rmse, 4),
      omega_rmse = round(omega_rmse, 4),
      AUC = AUC,
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
    pre_sigma  <- computeAssociations(model)
    sigma_hat  <- ((pre_sigma[[1]]$support > posterior_proba) + (pre_sigma[[1]]$support < (1-posterior_proba)) > 0) * pre_sigma[[1]]$mean
    omega_hat  <- solve(sigma_hat)
    omega_rmse <- Metrics::rmse(omega_hat, params$Omega)
    if(is.null(AUC)) AUC = get_auc(params$Omega, model, omega_hat, model_type)
    
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
      AUC = AUC,
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


    # Omega measures
    Sigma <- getResidualCov(model)$cov
    omega_hat <- solve(Sigma + diag(0.005, nrow(Sigma)))
    omega_rmse <- Metrics::rmse(omega_hat, params$Omega)
    if(is.null(AUC)) AUC = get_auc(params$Omega, model, omega_hat, model_type)

    # Fit measures
    fit_rmse   <- rmse(predict.gllvm(model), params$Y)

    # Aggregated results
    res <- c(
      criterion = NA,
      B_rmse = round(B_rmse, 4),
      omega_rmse = round(omega_rmse, 4),
      AUC = AUC,
      fit_rmse = round(fit_rmse, 4),
      roc_metrics(params$Omega, omega_hat)
    )
  }
  res
}