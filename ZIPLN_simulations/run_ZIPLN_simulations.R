############################ Loading useful libraries ##########################
source("ZIPLN_simulations.R")
set.seed(2)
distribution = "ZIP"

# Creating the saving repositories if needed
for (d in c("RES/ZIP_simulations",
            "RES/NB_simulations")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

############################ Reference zero-inflation values ###################
##### SITES #####
zi_mode_values_sites_ref = c(0.1, 0.25, 0.5, 0.75)
proba_mode_zi_values_sites_ref = rep(0.25, 4)

##### SPECIES #####
zi_mode_values_species_ref = c(0.1, 0.25, 0.5, 0.75)
proba_mode_zi_values_species_ref = rep(0.25, 4)


##### COVAR #####
block_values_ref_1 <- matrix(unlist(lapply(unlist(lapply(0.1 * 0.05 * seq(from = 2, to = 16), f <- function(x){max(x, 0.1)})), logit)), nrow = 3, byrow = TRUE)
block_values_ref_2 <- matrix(unlist(lapply(0.5 *0.05 * seq(from = 2, to = 16), logit)), nrow = 3, byrow = TRUE)
block_values_ref_3 <- matrix(unlist(lapply(0.05 * seq(from = 2, to = 16), logit)), nrow = 3, byrow = TRUE)
block_values_ref_4 <- matrix(unlist(lapply(1.1 * 0.05 * seq(from = 2, to = 16), logit)), nrow = 3, byrow = TRUE)

row_clusters_proba_ref <- rep(1/3, 3)
col_clusters_proba_ref <- rep(0.2, 5)

############################ Simulations parameters ############################
n_simu               <- 10
n_list               <- c(50)
p_list               <- c(25)
omega_structure_list <- c("erdos_renyi")

zi_type_list                <- c("species")
zi_mode_values_sites_list   <- NULL
proba_mode_zi_sites_list    <- NULL
zi_mode_values_species_list = list(0.1 * zi_mode_values_species_ref,
                                   0.5 * zi_mode_values_species_ref,
                                   zi_mode_values_species_ref,
                                   1.1 * zi_mode_values_species_ref
                                   )
proba_mode_zi_species_list = list(proba_mode_zi_values_species_ref,
                                  proba_mode_zi_values_species_ref,
                                  proba_mode_zi_values_species_ref,
                                  proba_mode_zi_values_species_ref
                                  )
block_values_list           <- NULL
row_clusters_proba_list     <- NULL
col_clusters_proba_list     <- NULL
################################################################################

res <- grid_ZIPLN_simulation(n_simu, n_list, p_list, omega_structure_list,
                             sigma_sparse = FALSE,
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

write.csv(res, 
          "RES/ZIP_simulations/res_1.csv")
