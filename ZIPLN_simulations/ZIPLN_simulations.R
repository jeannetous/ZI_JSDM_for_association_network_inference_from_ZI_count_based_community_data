library(PLNmodels)
library(tidyr)
library(dplyr)
library(Hmsc)
library(gllvm)
source("utils.R")
source("ZIPLN_measures.R")
source("ZIPLN_simulate_data.R")
source("negative_binomial_simulate_data.R")

#' @description simulates data under the ZIPLN model, runs different models and
#' outputs the performance measures for each one
#' @param simu rank of the simulation run ()
#' @param simu_params simulation parameters to simulate the data
#' @param PLN_formula formula to use to run PLN-network, the variables are named
#' V1, V2... so the formula should look like "Abundance ~ 0 + V1..."
#' @param ZIPLN_formula formula to use to run ZIPLN-network. If 
#' simu_params$zi_type = "covar", the ZI variables are named VZI1, VZI2...so the
#' formula should look like "Abundance ~ 0 + V1... | 0 + VZI1" If
#' simu_params$zi_type = "sites" or "species", the correct "zi" parameter is
#' used in ZIPLNnetwork
#' @param PLN_formula_ZIvar formula to use to run PLN-network including ZI
#' variables in the abundance. Useful only if zi_type = covar. Should look like
#' "Abundance ~ 0 + V1 +... + VZI1 + ..."
#' @param hmsc_params parameters for hmsc optimization
one_ZIPLN_simulation <- function(simu = 1, zi_config, simu_params,
                                 PLN_formula, ZIPLN_formula, 
                                 PLN_formula_ZIvar = NA,
                                 hmsc_params = NULL,
                                 distribution = c("ZIP", "negative.binomial")){
  distribution <- match.arg(distribution)
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
  ########################## Running hmsc  #####################################
  # Data formatting for hmsc
  Y_hurdle <- params$Y ; Y_hurdle[Y_hurdle == 0] = NA
  params$Y_hurdle <- Y_hurdle
  XData = data.frame("V1" = X[,"V1"])
  params$XData <- XData

  # hmsc parameters
  if(is.null(hmsc_params)){
    studyDesign = data.frame(sample = as.factor(1:simu_params$n))
    thin = 10
    hmsc_params <- list(studyDesign = studyDesign,
                        rL = HmscRandomLevel(units = studyDesign$sample),
                        posterior_proba = 0.9, thin = thin, samples= 1000, #500,
                        transient = 500 * thin, nChains = 1
    )
  }

  # Setting regression formula
  if(simu_params$add_intercept){
    XFormula <- "~ 1 + V1"
  }else{XFormula<- "~ V1 - 1"}

  # Running hmsc optimisation
  t0 = Sys.time()
  m = Hmsc(Y = Y_hurdle, XData = XData, XFormula = as.formula(XFormula),
           studyDesign = hmsc_params$studyDesign,
           ranLevels = list(sample = hmsc_params$rL),
           distr = "lognormal poisson")
  m = sampleMcmc(m, thin = hmsc_params$thin, samples = hmsc_params$samples,
                 transient = hmsc_params$transient,
                 nChains = hmsc_params$nChains, verbose = 0)
  t_hmsc = Sys.time() - t0
  hmsc_measures <- get_measures(m, params, model_selection = NULL,
                                stability = NULL,
                                AUC = NULL, add_intercept = simu_params$add_intercept,
                                model_type = "hmsc")
  hmsc_measures[["time"]] = as.numeric(t_hmsc, units = "secs")

  ############### Running hmsc with ZI covar, if applicable ####################
  if(!is.na(PLN_formula_ZIvar)){
    XData_ZIvar = data.frame("V1" = X[,"V1"], "VZI1" = X[,"VZI1"])
    params$XData_ZIvar <- XData_ZIvar
    t0 = Sys.time()
    m_ZIvar = Hmsc(Y = Y_hurdle, XData = XData_ZIvar, XFormula = ~ V1 + VZI1,
                   studyDesign = hmsc_params$studyDesign,
                   ranLevels = list(sample = hmsc_params$rL),
                   distr = "lognormal poisson")
    m_ZIvar = sampleMcmc(m, thin = hmsc_params$thin, samples = hmsc_params$samples,
                         transient = hmsc_params$transient,
                         nChains = hmsc_params$nChains, verbose = 0)
    t_hmsc_ZIvar = Sys.time() - t0
    hmsc_ZIvar_measures <- get_measures(m_ZIvar, params, model_selection = NULL,
                                        stability = NULL,
                                        add_intercept = simu_params$add_intercept,
                                        AUC = NULL, model_type = "hmsc")

    hmsc_ZIvar_measures[["time"]] = as.numeric(t_hmsc_ZIvar, units = "secs")
  }else{
    hmsc_ZIvar_measures <- NULL
  }
  
  ########################## Running gllvm #####################################
  if(simu_params$add_intercept){
    gllvm_formula <- "~ Intercept + V1"
  }else{gllvm_formula <- "~ V1"}
  t0 = Sys.time()
  gllvm_model <- gllvm(params$Y, params$X, family = distribution,
                       formula = as.formula(gllvm_formula))
  t_gllvm = Sys.time() - t0
  gllvm_measures <- get_measures(gllvm_model, params, model_selection = NULL,
                                stability = NULL,
                                add_intercept = simu_params$add_intercept,
                                AUC = NULL, model_type = "gllvm")
  gllvm_measures[["time"]] = as.numeric(t_gllvm, units = "secs")
  ############### Running glllvm with ZI covar, if applicable ##################
  if(!is.na(PLN_formula_ZIvar)){
    if(simu_params$add_intercept){
      gllvm_ZIvar_formula <- "~ Intercept + V1 + VZI1"
    }else{gllvm_ZIvar_formula <- "~ V1 + VZI1"}
    t0 = Sys.time()
    gllvm_ZIvar_model <- gllvm(params$Y, X, family = "ZIP", formula = as.formula(gllvm_ZIvar_formula))
    t_gllvm = Sys.time() - t0
    gllvm_ZIvar_measures <- get_measures(gllvm_ZIvar_model, params,
                                         model_selection = NULL, stability = NULL,
                                         add_intercept = simu_params$add_intercept,
                                         AUC = NULL, model_type = "gllvm")
    gllvm_ZIvar_measures[["time"]] = as.numeric(t_gllvm, units = "secs")
  }else{
    gllvm_ZIvar_measures <- NULL
  }
  
  ########################## Running PLN #######################################
  t0 = Sys.time()
  myPLN <- PLNnetwork(as.formula(PLN_formula), simu_data,
                      control = PLNnetwork_param(penalize_diagonal = FALSE,
                                                 min_ratio = 0.01,
                                                 n_penalties = 20,
                                                 trace = 0))
  t_PLN = Sys.time() - t0
  myPLN0 <- PLN(as.formula(PLN_formula), simu_data)
  PLN_StARS_1_measures <- get_measures(myPLN, params, model_selection = "StARS",
                                       unpenalised_model = myPLN0, stability = 0.8)
  t_PLN_StARS = Sys.time() - t0
  PLN_StARS_2_measures <- get_measures(myPLN, params, model_selection = "StARS",
                                       stability = 0.9,
                                       add_intercept = simu_params$add_intercept,
                                       AUC = PLN_StARS_1_measures[["AUC"]],
                                       AUC.alt = PLN_StARS_1_measures[["AUC.alt"]],
                                       Sigma_AUC = PLN_StARS_1_measures[["Sigma_AUC"]])
  
  PLN_BIC_measures <- get_measures(myPLN, params, model_selection = "BIC",
                                   add_intercept = simu_params$add_intercept,
                                   AUC = PLN_StARS_1_measures[["AUC"]],
                                   AUC.alt = PLN_StARS_1_measures[["AUC.alt"]],
                                   Sigma_AUC = PLN_StARS_1_measures[["Sigma_AUC"]])

  PLN_StARS_1_measures[["time"]] = as.numeric(t_PLN_StARS, units = "secs")
  PLN_StARS_2_measures[["time"]] = as.numeric(t_PLN_StARS, units = "secs")
  PLN_BIC_measures[["time"]] = as.numeric(t_PLN, units = "secs")
  
  ############### Running PLN with ZI covar, if applicable #####################
  if(!is.na(PLN_formula_ZIvar)){
    t0 = Sys.time()
    myPLN_ZIvar <- PLNnetwork(as.formula(PLN_formula_ZIvar), simu_data,
                              control = PLNnetwork_param(penalize_diagonal = FALSE,
                                                         min_ratio = 0.01,
                                                         n_penalties = 20))
    t_PLN_ZIvar = Sys.time() - t0
    myPLN_ZIvar0 <- PLN(as.formula(PLN_formula_ZIvar), simu_data)
    PLN_ZIvar_StARS_1_measures <- get_measures(myPLN_ZIvar, params, model_selection = "StARS",
                                               stability = 0.8, unpenalised_model = myPLN_ZIvar0,
                                               add_intercept = simu_params$add_intercept)
    t_PLN_ZIvar_StARS = Sys.time() - t0
    PLN_ZIvar_StARS_2_measures <- get_measures(myPLN_ZIvar, params,
                                               model_selection = "StARS",
                                               stability = 0.9,
                                               AUC = PLN_ZIvar_StARS_1_measures[["AUC"]],
                                               AUC.alt = PLN_ZIvar_StARS_1_measures[["AUC.alt"]],
                                               Sigma_AUC = PLN_ZIvar_StARS_1_measures[["Sigma_AUC"]])
    PLN_ZIvar_BIC_measures <- get_measures(myPLN_ZIvar, params, model_selection = "BIC",
                                           AUC = PLN_ZIvar_StARS_1_measures[["AUC"]],
                                           AUC.alt = PLN_ZIvar_StARS_1_measures[["AUC.alt"]],
                                           Sigma_AUC = PLN_ZIvar_StARS_1_measures[["Sigma_AUC"]])
    PLN_ZIvar_StARS_1_measures[["time"]] = as.numeric(t_PLN_ZIvar_StARS, units = "secs")
    PLN_ZIvar_StARS_2_measures[["time"]] = as.numeric(t_PLN_ZIvar_StARS, units = "secs")
    PLN_ZIvar_BIC_measures[["time"]] = as.numeric(t_PLN_ZIvar, units = "secs")
  }else{
    PLN_ZIvar_StARS_1_measures <- NULL ; PLN_ZIvar_StARS_2_measures <- NULL
    PLN_ZIvar_BIC_measures <- NULL
  }


  ######################### Running ZIPLN ######################################
  if(simu_params$zi_type == "covar"){zi <- NULL
  }else{zi <- ifelse(simu_params$zi_type == "sites", "row", "col")}

  t0 = Sys.time()
  myZIPLN <- ZIPLNnetwork(as.formula(ZIPLN_formula), simu_data, zi = zi,
                          control = ZIPLNnetwork_param(penalize_diagonal = FALSE,
                                                       min_ratio = 0.01,
                                                       n_penalties = 20,
                                                       trace = 0))

  myZIPLN0 <- ZIPLN(as.formula(ZIPLN_formula), simu_data, zi = zi)


  t_ZIPLN = Sys.time() - t0
  ZIPLN_StARS_1_measures <- get_measures(myZIPLN, params,
                                         unpenalised_model = myZIPLN0,
                                         model_selection = "StARS",
                                         stability = 0.8, model_type = "ZIPLN")
  t_ZIPLN_StARS = Sys.time() - t0
  ZIPLN_StARS_2_measures <- get_measures(myZIPLN, params, model_selection = "StARS",
                                         stability = 0.9,
                                         AUC = ZIPLN_StARS_1_measures[["AUC"]],
                                         AUC.alt = ZIPLN_StARS_1_measures[["AUC.alt"]],
                                         Sigma_AUC = ZIPLN_StARS_1_measures[["Sigma_AUC"]],
                                         model_type = "ZIPLN")
  ZIPLN_BIC_measures <- get_measures(myZIPLN, params, model_selection = "BIC",
                                     AUC = ZIPLN_StARS_1_measures[["AUC"]],
                                     AUC.alt = ZIPLN_StARS_1_measures[["AUC.alt"]],
                                     Sigma_AUC = ZIPLN_StARS_1_measures[["Sigma_AUC"]],
                                     model_type = "ZIPLN")
  ZIPLN_StARS_1_measures[["time"]] = as.numeric(t_ZIPLN_StARS, units = "secs")
  ZIPLN_StARS_2_measures[["time"]] = as.numeric(t_ZIPLN_StARS, units = "secs")
  ZIPLN_BIC_measures[["time"]] = as.numeric(t_ZIPLN, units = "secs")
  
  ################# Merging all the measures in one data frame #################
  measure_rows <- list(c(method = "PLN", PLN_StARS_1_measures),
                       c(method = "PLN", PLN_StARS_2_measures),
                       c(method = "PLN", PLN_BIC_measures),
                       c(method = "ZIPLN", ZIPLN_StARS_1_measures),
                       c(method = "ZIPLN", ZIPLN_StARS_2_measures),
                       c(method = "ZIPLN", ZIPLN_BIC_measures),
                       c(method = "hmsc", hmsc_measures),
                       c(method = "gllvm", gllvm_measures)
                       )

  if(!is.na(PLN_formula_ZIvar)){
    measure_rows <- c(measure_rows,
                      list(c(method = "PLN_ZIvar", PLN_ZIvar_StARS_1_measures),
                           c(method = "PLN_ZIvar", PLN_ZIvar_StARS_2_measures),
                           c(method = "PLN_ZIvar", PLN_ZIvar_BIC_measures),
                           c(method = "hmsc_ZIvar", hmsc_ZIvar_measures),
                           c(method = "gllvm_ZIvar", gllvm_ZIvar_measures)))
  }
  

  res <- as.data.frame(cbind(simu = simu, n = simu_params$n, p = simu_params$p,
                             omega_structure = simu_params$omega_structure,
                             distribution = distribution,
                             zi_type = simu_params$zi_type,
                             zi_config = zi_config,
                             do.call(rbind, measure_rows)))
  return(res)
}


#' @description runs multiple simulations with the same parameters
#' @param n_simu number of simulations to run
#' #' @param PLN_formula formula to use to run PLN-network, the variables are named
#' V1, V2... so the formula should look like "Abundance ~ 0 + V1..."
#' @param ZIPLN_formula formula to use to run ZIPLN-network. If simu_params$zi_type = "covar",
#' the ZI variables are named VZI1, VZI2...so the formula should look like "Abundance ~ 0 + V1... | 0 + VZI1"
#' If simu_params$zi_type = "sites" or "species", the correct "zi" parameter is used in ZIPLNnetwork
#' @param PLN_formula_ZIvar formula to use to run PLN-network including ZI variables in the abundance. Useful
#' only if zi_type = covar. Should look like "Abundance ~ 0 + V1 +... + VZI1 + ..."
#' @param mc.cores number of cores to run the simulations on in parallel
multiple_ZIPLN_simulations <- function(n_simu, zi_config, simu_params,
                                       PLN_formula, ZIPLN_formula,
                                       PLN_formula_ZIvar = NA,
                                       hmsc_params = NULL, distribution,
                                       mc.cores = max(1, parallel::detectCores() - 2)){
  cat("Settings: (n, p, omega structure, zi type) = (",simu_params$n, simu_params$p,
      simu_params$omega_structure, simu_params$zi_type, ")\n")

  multiple_res <- parallel::mclapply(1:n_simu,
                                     one_ZIPLN_simulation,
                                     zi_config = zi_config,
                                     simu_params = simu_params,
                                     PLN_formula = PLN_formula,
                                     ZIPLN_formula = ZIPLN_formula,
                                     PLN_formula_ZIvar = PLN_formula_ZIvar,
                                     hmsc_params = hmsc_params,
                                     distribution = distribution,
                                     mc.cores = mc.cores)
  res <- do.call(rbind, multiple_res)
  res
}

#' @description runs multiple simulations for each set of parameters in a grid
#' @param n_simu number of simulations to run for each parametrization
#' @param n_list number of n (number of sites) values to go through
#' @param p_list number of p (number of species) values to go through
#' @param omega_structure_list list of omega_structure values to go through
#' @param sigma_sparse  boolean, indicates whether the simulation should be
#' based on a sparse variance matrix - default is FALSE, corresponding to a
#' sparse precision matrix
#' @param zi_type_list list of zi_type values to go through
#' @param zi_mode_values_sites_list, list of zi_mode_values values to go through for the sites
#' @param proba_mode_zi_sites_list, list of proba_mode_zi values to go through for the sites
#' @param zi_mode_values_species_list, list of zi_mode_values values to go through for the species
#' @param proba_mode_zi_species_list, list of proba_mode_zi values to go through for the species
#' @param block_values_list list of matrices for block_values for when zi_type = covar
#' @param row_clusters_proba_list list of lists of probabilities for row clusters, used only if row_clusters=NULL,
#' default is equiprobable distributions
#' @param col_clusters_proba_list list of lists of probabilities for column clusters, used only if col_clusters=NULL,
#' default is equiprobable distributions
#' @param min_X minimum value for X, either one single value for X, or a list of
#' length d for each dimension, fixed along the grid
#' @param max_X maximum value for X, either one single value for X, or a list of
#' length d for each dimension, fixed along the grid
#' @param mean_B mean of the Gaussian distributions B values are drawn from, can be a list with one mean per B dimension
#' @param sd_B standard deviation of the Gaussian distributions B values are drawn from, can be a list with one sd per B dimension
#' @param XB_max max value for exp(XB) mean of the Poisson distribution that Y follows
#' @param min_X0 minimum value for X0, either one single value for X0, or a
#' list of length d for each dimension, applied only if zi_covar_cluster = FALSE
#' (otherwise X0 is discrete), fixed along the grid
#' @param max_X0 maximum value for X0, either one single value for X0, or a
#' list of length d for each dimension applied only if zi_covar_cluster = FALSE
#' (otherwise X0 is discrete), fixed along the grid
#' @param max_X0B0 maximum value for the mean of each column of X0 %*% B0
#' applied only if zi_covar_cluster = FALSE (otherwise X0 is discrete), fixed
#' along the grid
#' @param hmsc_params parameters for hmsc optimization
#' @param mc.cores number of cores to run the simulations on in parallel
grid_ZIPLN_simulation <- function(n_simu, n_list, p_list, add_intercept,
                                  omega_structure_list, sigma_sparse = FALSE,
                                  zi_type_list, zi_mode_values_sites_list,
                                  proba_mode_zi_sites_list,
                                  zi_mode_values_species_list,
                                  proba_mode_zi_species_list, block_values_list,
                                  row_clusters_proba_list, col_clusters_proba_list,
                                  min_X = 0,  max_X = 1, mean_B  = 2, sd_B = 1,
                                  XB_max = 70, min_X0 = 0, max_X0 = 10,
                                  max_X0B0 = -0.2, d = 1, hmsc_params = NULL,
                                  distribution = distribution,
                                  mc.cores = max(1, parallel::detectCores() - 2)) {
  settings <- expand.grid(n = n_list,
                          p = p_list,
                          omega_structure = omega_structure_list,
                          zi_type = zi_type_list,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  if(!is.null(zi_mode_values_sites_list)){
    sites_zi_tuples <- tibble(zi_mode_values  = zi_mode_values_sites_list,
                              proba_mode_zi   = proba_mode_zi_sites_list,
                              zi_config = paste0("sites_", 1:length(zi_mode_values_sites_list)))
  }
  if(!is.null(zi_mode_values_species_list)){
    species_zi_tuples <- tibble(zi_mode_values  = zi_mode_values_species_list,
                                proba_mode_zi   = proba_mode_zi_species_list,
                                zi_config = paste0("species_", 1:length(zi_mode_values_species_list)))
  }
  
  
  settings <- settings %>%
    rowwise() %>%
    do({
      row <- .
      if (row$zi_type == "covar") {
        tibble(n = row$n, p = row$p, omega_structure = row$omega_structure,
               zi_type = row$zi_type, zi_mode_values = NA,
               proba_mode_zi = NA)
      } else if (row$zi_type == "sites") {
        cbind(row[1:4], sites_zi_tuples)
      } else if (row$zi_type == "species") {
        cbind(row[1:4], species_zi_tuples)
      }
    }) %>%
    ungroup()
  
  if(! is.null(block_values_list)){
    block_values_rows <- do.call(rbind, lapply(seq_along(block_values_list), function(i) {
      df_tmp <- settings %>% filter(zi_type == "covar")
      df_tmp$block_values <- rep(list(block_values_list[[i]]), nrow(df_tmp))
      df_tmp$row_clusters_proba <- rep(list(row_clusters_proba_list[[i]]), nrow(df_tmp))
      df_tmp$col_clusters_proba <- rep(list(col_clusters_proba_list[[i]]), nrow(df_tmp))
      df_tmp$zi_config <- rep(paste0("covar_", i), nrow(df_tmp))
      return(df_tmp)
    }))
    
    settings <- settings %>% filter(zi_type != "covar")
    settings <- settings %>% mutate(block_values = NA, row_clusters_proba = NA, col_clusters_proba = NA)
    settings <- rbind(settings, block_values_rows)
  }else{settings$block_values <- NA
  settings$row_clusters_proba <- NA ; settings$col_clusters_proba <- NA}
  
  settings$PLN_formula_ZIvar <- NA
  if(add_intercept){
    settings$PLN_formula <- "Abundance ~ 1 + V1"
    settings$ZIPLN_formula <- "Abundance ~ 1 + V1"
    settings[settings$zi_type == "covar",]$ZIPLN_formula <- "Abundance ~ 1 + V1 | 0 + VZI1"
    settings[settings$zi_type == "covar",]$PLN_formula_ZIvar <- "Abundance ~ 1 + V1 + VZI1"
  }else{
    settings$PLN_formula <- "Abundance ~ 0 + V1"
    settings$ZIPLN_formula <- "Abundance ~ 0 + V1"
    settings[settings$zi_type == "covar",]$ZIPLN_formula <- "Abundance ~ 0 + V1 | 0 + VZI1"
    settings[settings$zi_type == "covar",]$PLN_formula_ZIvar <- "Abundance ~ 0 + V1 + VZI1"
  }
  settings$add_intercept <- add_intercept
  settings$n_simu <- n_simu
  final_res <- purrr::pmap(settings, f <- function(n, p, omega_structure, zi_type,
                                                   zi_mode_values,
                                                   proba_mode_zi, zi_config,
                                                   block_values,
                                                   row_clusters_proba,
                                                   col_clusters_proba,
                                                   PLN_formula,
                                                   ZIPLN_formula, PLN_formula_ZIvar,
                                                   add_intercept,
                                                   n_simu){
    
    simu_params = list(n = n,
                       p = p,
                       d = d,
                       add_intercept = add_intercept,
                       omega_structure = omega_structure,
                       sigma_sparse    = sigma_sparse,
                       zi_type = zi_type,
                       zi_covar_cluster = TRUE,
                       min_X = min_X, max_X = max_X,
                       mean_B  = mean_B, sd_B = sd_B, XB_max = XB_max,
                       zi_mode_values = zi_mode_values,
                       proba_mode_zi = proba_mode_zi,
                       block_values = block_values,
                       row_clusters = NULL,
                       col_clusters = NULL,
                       row_clusters_proba = row_clusters_proba,
                       col_clusters_proba = col_clusters_proba,
                       X0 = NULL, B0 = NULL,
                       min_X0 = min_X0, max_X0 = max_X0,
                       max_X0B0 = max_X0B0)
    
    multiple_ZIPLN_simulations(
      n_simu = n_simu, zi_config = zi_config, simu_params,
      PLN_formula = PLN_formula,
      ZIPLN_formula = ZIPLN_formula,
      PLN_formula_ZIvar = PLN_formula_ZIvar,
      hmsc_params = hmsc_params,
      distribution = distribution
    )})
  final_res <- final_res[sapply(final_res, function(df) ncol(df) > 1)]
  final_res <- do.call(rbind, final_res) %>% as_tibble()
  final_res
}
