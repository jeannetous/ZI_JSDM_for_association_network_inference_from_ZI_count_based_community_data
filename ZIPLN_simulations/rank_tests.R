source("ZIPLN_simulations.R")

# Creating the saving repository if needed
if (!dir.exists("Rank_tests")) dir.create("Rank_tests", recursive = TRUE)

########################### SIMULATION PARAMETERS ##############################
distribution = "ZIP"
zi_type = "species"
n_simu = 1
rank_list = 2:3
zi_mode_values_species_ref = c(0.1, 0.25, 0.5, 0.75)
proba_mode_zi_values_species_ref = rep(0.25, 4)

simu_params = list(n = 100,
                   p = 25,
                   d = 1,
                   sigma_sparse = TRUE,
                   add_intercept = TRUE,
                   omega_structure = "erdos_renyi",
                   zi_type = zi_type,
                   zi_covar_cluster = TRUE,
                   min_X = 0, max_X = 1, mean_B  = 2, sd_B = 1, XB_max = 70,
                   zi_mode_values = zi_mode_values_species_ref,
                   proba_mode_zi = proba_mode_zi_values_species_ref, 
                   block_values = NULL, 
                   row_clusters = NULL,
                   col_clusters = NULL,
                   row_clusters_proba = NULL, 
                   col_clusters_proba = NULL, 
                   X0 = NULL, B0 = NULL,
                   min_X0 = 0, max_X0 = 10,
                   max_X0B0 = 0.2)


########################## gllvm ###############################################
res <- list()
for(simu in 1:n_simu){
  cat(paste0("SIMU #", as.character(simu)))
  Y <- matrix(rep(0, simu_params$n, simu_params$p), nrow = simu_params$n)
  while( (TRUE %in% (colSums(Y) == 0)) | (TRUE %in% (rowSums(Y) == 0)) ){
    params <- do.call(generate_all_ZIPLN_parameters, simu_params)
    if(distribution == "ZIP"){
      Y <- simulate_ZIPLN_data(params)
    }else{
      params$r <- 10
      Y <- simulate_negative_binomial_data(params)}
  }
  if(!is.null(params$zi_params)){
    X <- data.frame(params$X, as.factor(params$zi_params$X0))
    colnames(X)[[length(colnames(X))]] <- "VZI1"
  }else{X <- params$X}
  simu_data <- prepare_data(Y, X)
  params$Y <- simu_data$Abundance

  if(simu_params$add_intercept){
    gllvm_formula <- "~ Intercept + V1"
  }else{gllvm_formula <- "~ V1"}

  for(i in rank_list){
    print(paste0("num.lv = ", as.character(i)))
    t0 = Sys.time()
    gllvm_model <- gllvm(params$Y, params$X, family = distribution,
                         num.lv = i,
                         formula = as.formula(gllvm_formula))
    t_gllvm = Sys.time() - t0
    gllvm_measures <- get_measures(gllvm_model, params, model_selection = NULL,
                                   stability = NULL,
                                   add_intercept = simu_params$add_intercept,
                                   AUC = NULL, model_type = "gllvm")
    gllvm_measures[["num.lv"]]      = gllvm_model$num.lv
    gllvm_measures[["convergence"]] = gllvm_model$convergence
    gllvm_measures[["time"]]        = as.numeric(t_gllvm, units = "secs")
    
    res <- rbind(res, gllvm_measures)
  }
}

write.csv(res, "Rank_tests/gllvm_rank.csv")

########################## hmsc ################################################
# hmsc parameters
studyDesign = data.frame(sample = as.factor(1:simu_params$n))
thin = 10
hmsc_params <- list(studyDesign = studyDesign,
                    rL = HmscRandomLevel(units = studyDesign$sample),
                    posterior_proba = 0.9, thin = thin, samples= 1000, 
                    transient = 500 * thin, nChains = 2)

res <- list()

for(simu in 1:n_simu){
  cat(paste0("SIMU #", as.character(simu)))
  Y <- matrix(rep(0, simu_params$n, simu_params$p), nrow = simu_params$n)
  while( (TRUE %in% (colSums(Y) == 0)) | (TRUE %in% (rowSums(Y) == 0)) ){
    params <- do.call(generate_all_ZIPLN_parameters, simu_params)
    if(distribution == "ZIP"){
      Y <- simulate_ZIPLN_data(params)
    }else{
      params$r <- 10
      Y <- simulate_negative_binomial_data(params)}

  }
  if(!is.null(params$zi_params)){
    X <- data.frame(params$X, as.factor(params$zi_params$X0))
    colnames(X)[[length(colnames(X))]] <- "VZI1"
  }else{X <- params$X}
  simu_data <- prepare_data(Y, X)
  params$Y <- simu_data$Abundance

  # Data formatting for hmsc
  Y_hurdle <- params$Y ; Y_hurdle[Y_hurdle == 0] = NA
  params$Y_hurdle <- Y_hurdle
  XData = data.frame("V1" = X[,"V1"])
  params$XData <- XData

  # Setting regression formula
  if(simu_params$add_intercept){
    XFormula <- "~ 1 + V1"
  }else{XFormula<- "~ V1 - 1"}

  for(i in rank_list){
    print(paste0("num.lv = ", as.character(i)))
    rL <- HmscRandomLevel(units = studyDesign$sample)
    rL <- setPriors(rL, nfMin = i, nfMax = i)

    m <- Hmsc(Y = Y_hurdle, XData = XData, XFormula = as.formula(XFormula),
              studyDesign = hmsc_params$studyDesign,
              ranLevels = list(sample = rL), distr = "lognormal poisson")

    t0 = Sys.time()
    m = sampleMcmc(m, thin = hmsc_params$thin, samples = hmsc_params$samples,
                   transient = hmsc_params$transient,
                   nChains = hmsc_params$nChains, verbose = 0)
    t_hmsc = Sys.time() - t0
    mpost <- convertToCodaObject(m)
    convergence_check <- (abs(median(gelman.diag(mpost$Beta, multivariate = FALSE)$psrf[,1]) - 1) < 0.1
                          & abs(median(gelman.diag(mpost$V, multivariate = FALSE)$psrf[,1]) - 1) < 0.1)
    hmsc_measures <- get_measures(m, params, model_selection = NULL,
                                  stability = NULL, 
                                  add_intercept = simu_params$add_intercept,
                                  AUC = NULL, model_type = "hmsc")
    hmsc_measures[["num.lv"]] = dim(m$postList[[1]][[1]]$Lambda[[1]])[[1]]
    hmsc_measures[["convergence"]] = convergence_check
    hmsc_measures[["time"]] = as.numeric(t_hmsc, units = "secs")
    res <- rbind(res, hmsc_measures)
  }
}

write.csv(res, "Rank_tests/hmsc_rank.csv")
