############################ Loading useful libraries ##########################
source("ZIPLN_simulations.R")
set.seed(2)
distribution = "ZIP"
############################ Reference ZI values from real data - ZIPLN ########
##### SITES #####
zi_mode_values_sites_ref = c(0.2, 0.4, 0.6)
proba_mode_zi_values_sites_ref = c(0.3, 0.3, 0.4)

##### SPECIES #####
zi_mode_values_species_ref = c(0.05, 0.3, 0.6, 0.9)
proba_mode_zi_values_species_ref = c(0.4, 0.15, 0.15, 0.3)

##### SITES X SPECIES #####
block_values_ref0 <- matrix(-20, nrow = 5, ncol = 12)
block_values_ref <- matrix(c(-0.8, 1.5, 3.2, -0.8, -2,
                             -2.2, 7, -1, -2.7, 7,
                             -0.7, 7, -1, -0.6, -2,
                             -0.7, 7, 7, -3.6, -2,
                              7, 7, 7, -1.7, -2,
                              7, 7, 0.8, -0.4, -2,
                              -2, -4.8, -4.4, -1.9, -2,
                              7, -2, 7, -1, -2,
                              7, 7, 7, 0.6, -2,
                              7, 7, 7, -0.1, -2,
                              -0.1, 7, 7, 0, -2,
                              7, 7, 7, 0.2, -2 ), nrow = 5)
block_values_ref_1 <- block_values_ref ; block_values_ref_1[block_values_ref_1 > 0] <- 0.1 * block_values_ref_1[block_values_ref_1 > 0]
block_values_ref_2 <- block_values_ref ; block_values_ref_2[block_values_ref_2 > 0] <- 0.5 * block_values_ref_2[block_values_ref_2 > 0]
block_values_ref_3 <- block_values_ref
block_values_ref_4 <- block_values_ref ; block_values_ref_4[block_values_ref_4 > 0] <- 1.1 * block_values_ref_4[block_values_ref_4 > 0]
block_values_ref_4[block_values_ref_4 < 0] <- 0.5 * block_values_ref_4[block_values_ref_4 < 0]

row_clusters_proba_ref <- c(0.125, 0.125, 0.125, 0.125, 0.5)
col_clusters_proba_ref <- rep(0.0833, 12)

############################ Simulations parameters ############################
n_simu =  1
n_list = c(50)
p_list = c(25)
omega_structure_list = c("community")
zi_type_list =  c("species")


zi_mode_values_sites_list = NULL
proba_mode_zi_sites_list = NULL

zi_mode_values_species_list = list(zi_mode_values_species_ref)

proba_mode_zi_species_list = list(proba_mode_zi_values_species_ref)



block_values_list = NULL
row_clusters_proba_list = NULL ; col_clusters_proba_list = NULL

############################ Species-dep ZI ####################################
############################ Running simulations ###############################
res <- grid_ZIPLN_simulation(n_simu, n_list, p_list, omega_structure_list,
                             add_intercept = TRUE,
                             zi_type_list,
                             zi_mode_values_sites_list = zi_mode_values_sites_list,
                             proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                             zi_mode_values_species_list = zi_mode_values_species_list,
                             proba_mode_zi_species_list = proba_mode_zi_species_list,
                             block_values_list = block_values_list,
                             row_clusters_proba_list = row_clusters_proba_list,
                             col_clusters_proba_list = col_clusters_proba_list,
                             min_X = 0,  max_X = 1,
                             mean_B  = 2, sd_B = 1, XB_max = 70,
                             min_X0 = 0, max_X0 = 10, max_X0B0 = -0.2,
                             distribution = distribution, d = 1,
                             mc.cores = max(1, parallel::detectCores() - 2))

############################ Saving the results and parameters #################
write.csv(res, "simulations_res_species.csv")



all_params <- list(n_simu = n_simu, n_list = n_list, p_list = p_list,
                   min_X = 0,  max_X = 1,
                   mean_B  = 2, sd_B = 1, XB_max = 70,
                   omega_structure_list = omega_structure_list,
                   zi_type_list = zi_type_list,
                   zi_mode_values_sites_list = zi_mode_values_sites_list,
                   proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                   zi_mode_values_species_list = zi_mode_values_species_list,
                   proba_mode_zi_species_list = proba_mode_zi_species_list,
                   distribution = distribution)

writeLines(capture.output(str(all_params)), "simulations_res_species.txt")

############################ Sites-dep ZI ######################################
zi_type_list =  c("sites")
zi_mode_values_sites_list = list(zi_mode_values_sites_ref)

proba_mode_zi_sites_list = list(proba_mode_zi_values_sites_ref)
zi_mode_values_species_list = NULL ; proba_mode_zi_species_list = NULL
############################ Running simulations ###############################
res <- grid_ZIPLN_simulation(n_simu, n_list, p_list, omega_structure_list,
                             add_intercept = TRUE,
                             zi_type_list,
                             zi_mode_values_sites_list = zi_mode_values_sites_list,
                             proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                             zi_mode_values_species_list = zi_mode_values_species_list,
                             proba_mode_zi_species_list = proba_mode_zi_species_list,
                             block_values_list = block_values_list,
                             row_clusters_proba_list = row_clusters_proba_list,
                             col_clusters_proba_list = col_clusters_proba_list,
                             min_X = 0,  max_X = 1,
                             mean_B  = 2, sd_B = 1, XB_max = 70,
                             min_X0 = 0, max_X0 = 10, max_X0B0 = -0.2,
                             distribution = distribution,
                             mc.cores = max(1, parallel::detectCores() - 2))

############################ Saving the results and parameters #################
write.csv(res, "simulations_res_sites.csv")



all_params <- list(n_simu = n_simu, n_list = n_list, p_list = p_list,
                   min_X = 0,  max_X = 1,
                   mean_B  = 2, sd_B = 1, XB_max = 70,
                   omega_structure_list = omega_structure_list,
                   zi_type_list = zi_type_list,
                   zi_mode_values_sites_list = zi_mode_values_sites_list,
                   proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                   zi_mode_values_species_list = zi_mode_values_species_list,
                   proba_mode_zi_species_list = proba_mode_zi_species_list,
                   distribution = distribution)

writeLines(capture.output(str(all_params)), "simulations_res_sites_parameters.txt")



############################ Covar-dep ZI ######################################
zi_type_list =  c("covar")
block_values_list = list(block_values_ref_1)

row_clusters_proba_list = list(row_clusters_proba_ref)
col_clusters_proba_list = list(col_clusters_proba_ref)

############################ Running simulations ###############################
res <- grid_ZIPLN_simulation(n_simu, n_list, p_list, omega_structure_list,
                             add_intercept = TRUE,
                             zi_type_list,
                             zi_mode_values_sites_list = zi_mode_values_sites_list,
                             proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                             zi_mode_values_species_list = zi_mode_values_species_list,
                             proba_mode_zi_species_list = proba_mode_zi_species_list,
                             block_values_list = block_values_list,
                             row_clusters_proba_list = row_clusters_proba_list,
                             col_clusters_proba_list = col_clusters_proba_list,
                             min_X = 0,  max_X = 1,
                             mean_B  = 2, sd_B = 1, XB_max = 70,
                             min_X0 = 0, max_X0 = 10, max_X0B0 = -0.2,
                             distribution = distribution,
                             mc.cores = max(1, parallel::detectCores() - 2))

############################ Saving the results and parameters #################
write.csv(res, "simulations_res_covar.csv")



all_params <- list(n_simu = n_simu, n_list = n_list, p_list = p_list,
                   min_X = 0,  max_X = 1,
                   mean_B  = 2, sd_B = 1, XB_max = 70,
                   omega_structure_list = omega_structure_list,
                   zi_type_list = zi_type_list,
                   zi_mode_values_sites_list = zi_mode_values_sites_list,
                   proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                   zi_mode_values_species_list = zi_mode_values_species_list,
                   proba_mode_zi_species_list = proba_mode_zi_species_list,
                   distribution = distribution)

writeLines(capture.output(str(all_params)), "simulations_res_covar_parameters.txt")




############################ SIMULATIONS WITH UNOBSERVED COVARIATE #############

res <- grid_ZIPLN_simulation(n_simu, n_list, p_list, omega_structure_list,
                             add_intercept = TRUE,
                             zi_type_list,
                             zi_mode_values_sites_list = zi_mode_values_sites_list,
                             proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                             zi_mode_values_species_list = zi_mode_values_species_list,
                             proba_mode_zi_species_list = proba_mode_zi_species_list,
                             block_values_list = block_values_list,
                             row_clusters_proba_list = row_clusters_proba_list,
                             col_clusters_proba_list = col_clusters_proba_list,
                             min_X = 0,  max_X = 1,
                             mean_B  = 2, sd_B = 1, XB_max = 70,
                             min_X0 = 0, max_X0 = 10, max_X0B0 = -0.2,
                             distribution = distribution, d = 2,
                             mc.cores = max(1, parallel::detectCores() - 2))

############################ Saving the results and parameters #################
write.csv(res, "simulations_res_species_unobserved_covariate.csv")



all_params <- list(n_simu = n_simu, n_list = n_list, p_list = p_list,
                   min_X = 0,  max_X = 1,
                   mean_B  = 2, sd_B = 1, XB_max = 70,
                   omega_structure_list = omega_structure_list,
                   zi_type_list = zi_type_list,
                   zi_mode_values_sites_list = zi_mode_values_sites_list,
                   proba_mode_zi_sites_list = proba_mode_zi_sites_list,
                   zi_mode_values_species_list = zi_mode_values_species_list,
                   proba_mode_zi_species_list = proba_mode_zi_species_list,
                   distribution = distribution)

writeLines(capture.output(str(all_params)), "simulations_res_species_unobserved_covariate_parameters.txt")

