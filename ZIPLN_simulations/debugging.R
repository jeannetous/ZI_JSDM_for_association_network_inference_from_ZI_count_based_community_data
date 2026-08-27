set.seed(1)
source("ZIPLN_simulations.R")
n = 50 ; p = 25
min_X = 0;   max_X = 10; SNR = 0.75
mean_B  = 2 ; sd_B = 1 ; XB_max = 70
omega_structure = "erdos_renyi"
zi_type = "species"
zi_mode_values = c(0.1, 0.25, 0.5, 0.75)
proba_mode_zi  = rep(0.25, 4)

simu_params = list(n = n,
                   p = p,
                   d = 1,
                   add_intercept = TRUE,
                   omega_structure = omega_structure,
                   zi_type = zi_type,
                   zi_covar_cluster = TRUE,
                   min_X = 0, max_X = 1, mean_B  = 2, sd_B = 1, XB_max = 70,
                   zi_mode_values = zi_mode_values,
                   proba_mode_zi = proba_mode_zi, 
                   block_values = NULL, 
                   row_clusters = NULL,
                   col_clusters = NULL,
                   row_clusters_proba = NULL, 
                   col_clusters_proba = NULL, 
                   X0 = NULL, B0 = NULL,
                   min_X0 = 0, max_X0 = 10,
                   max_X0B0 = 0.2)

PLN_formula <- "Abundance ~ 1 + V1"
ZIPLN_formula <- "Abundance ~ 1 + V1"
PLN_formula_ZIvar <- NA 
res <- one_ZIPLN_simulation(1, "covar", simu_params, PLN_formula,
                            ZIPLN_formula, PLN_formula_ZIvar,
                            distribution = "ZIP")
