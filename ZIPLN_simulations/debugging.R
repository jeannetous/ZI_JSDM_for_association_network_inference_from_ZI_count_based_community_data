set.seed(1)
source("ZIPLN_simulations.R")
n = 300 ; p = 20
min_X = 0;   max_X = 10; SNR = 0.75
mean_B  = 2 ; sd_B = 1 ; XB_max = 70
omega_structure = "erdos_renyi"

# zi_type = "sites"
# n_mode_zi_proba = 3
# zi_mode_values =  c(0.2, 0.4, 0.6)
# proba_mode_zi = c(0.3, 0.3, 0.4)

# zi_type = "species"
# n_mode_zi_proba = 4
# zi_mode_values = c(0.2, 0.3, 0.6, 0.95)
# Les deux lignes du dessous ne donnent pas lieu à de très bonnes AUC
# zi_mode_values = c(0.05, 0.3, 0.6, 0.9)
# proba_mode_zi = c(0.3, 0.2, 0.15, 0.35)
# zi_mode_values = c(0.05, 0.3, 0.6, 0.9)
# proba_mode_zi = c(0.4, 0.15, 0.15, 0.3)
# zi_mode_values = c(0.05, 0.1, 0.15, 0.2)
# zi_mode_values = c(0.6, 0.65, 0.7, 0.75)
# proba_mode_zi = c(0.1, 0.2, 0.3, 0.4)

# zi_type = "species"
# block_values_ref <- matrix(c(-0.8, 1.5, 3.2, -0.8, -2,
#                              -2.2, 7, -1, -2.7, 7,
#                              -0.7, 7, -1, -0.6, -2,
#                              -0.7, 7, 7, -3.6, -2,
#                              7, 7, 7, -1.7, -2,
#                              7, 7, 0.8, -0.4, -2,
#                              -2, -4.8, -4.4, -1.9, -2,
#                              7, -2, 7, -1, -2,
#                              7, 7, 7, 0.6, -2,
#                              7, 7, 7, -0.1, -2,
#                              -0.1, 7, 7, 0, -2,
#                              7, 7, 7, 0.2, -2 ), nrow = 5)
# block_values_ref_1 <- block_values_ref ; block_values_ref_1[block_values_ref_1 > 0] <- 0.1 * block_values_ref_1[block_values_ref_1 > 0]
# block_values_ref_2 <- block_values_ref ; block_values_ref_2[block_values_ref_2 > 0] <- 0.5 * block_values_ref_2[block_values_ref_2 > 0]
# block_values_ref_3 <- block_values_ref
# block_values_ref_4 <- block_values_ref ; block_values_ref_4[block_values_ref_4 > 0] <- 1.1 * block_values_ref_4[block_values_ref_4 > 0]
# block_values_ref_4[block_values_ref_4 < 0] <- 0.5 * block_values_ref_4[block_values_ref_4 < 0]
# 
# row_clusters_proba_ref <- c(0.125, 0.125, 0.125, 0.125, 0.5)
# col_clusters_proba_ref <- rep(0.0833, 12)

zi_type = "species"
n_mode_zi_proba = 4
zi_mode_values_species_ref = c(0.05, 0.3, 0.6, 0.9)
zi_mode_values = 1.1 * zi_mode_values_species_ref
proba_mode_zi = c(0.4, 0.15, 0.15, 0.3)

simu_params = list(n = n,
                   p = p,
                   d = 1,
                   add_intercept = TRUE,
                   omega_structure = omega_structure,
                   zi_type = zi_type,
                   zi_covar_cluster = TRUE,
                   min_X = 0, max_X = 1, mean_B  = 2, sd_B = 1, XB_max = 70,
                   zi_mode_values = zi_mode_values,#NULL, #
                   proba_mode_zi = proba_mode_zi, #NULL, #
                   block_values = NULL, #block_values_ref, #0.5 * # #
                   row_clusters = NULL,
                   col_clusters = NULL,
                   row_clusters_proba = NULL, #row_clusters_proba_ref,#
                   col_clusters_proba = NULL, #col_clusters_proba_ref,#
                   X0 = NULL, B0 = NULL,
                   min_X0 = 0, max_X0 = 10,
                   max_X0B0 = 0.2)

# res <- multiple_ZIPLN_simulations(3, "sites_1", simu_params, "Abundance ~ 0 + V1", "Abundance ~ 0 + V1")
# for(i in 1:10){
  # print(i)
  i = 1
  set.seed(i)
  PLN_formula <- "Abundance ~ 1 + V1"
  ZIPLN_formula <- "Abundance ~ 1 + V1"# "Abundance ~ 1 + V1 | 0 + VZI1" #
  PLN_formula_ZIvar <- NA # "Abundance ~ 1 + V1 + VZI1" #
  simu_params$d = 2
  res <- one_ZIPLN_simulation(1, "covar", simu_params, PLN_formula,
                              ZIPLN_formula, PLN_formula_ZIvar,
                              distribution = "ZIP") #| 0 +  VZI1
                              # "Abundance ~ 1 + V1 + VZI1")
  # res <- multiple_ZIPLN_simulations(30, "covar", simu_params, "Abundance ~ 1 + V1",
  #                                   "Abundance ~ 1 + V1", NA) # | 0 +  VZI1"
  # "Abundance ~ 1 + V1 + VZI1")
  # print(res)
# }
# res <- one_ZIPLN_simulation(1, "covar_2", simu_params, "Abundance ~ 1 + V1",
#                             "Abundance ~ 1 + V1 | VZI1",
#                             "Abundance ~ 1 + V1 + VZI1")
# res <- one_ZIPLN_simulation(1, "sites_1", simu_params, "Abundance ~ 0 + V1", "Abundance ~ 0 + V1")


  ################################################################################


  # simu_params = list(n = 300,
  #                    p = 20,
  #                    d = 1,
  #                    omega_structure = "erdos_renyi",
  #                    zi_type = "species",
  #                    zi_covar_cluster = TRUE,
  #                    min_X = 0, max_X = 1, mean_B  = 2, sd_B = 1, XB_max = 70,
  #                    zi_mode_values = 0.1 * zi_mode_values_species_ref,#NULL, #
  #                    proba_mode_zi = proba_mode_zi_values_species_ref, #NULL, #
  #                    block_values = NULL, #block_values_ref, #0.5 *
  #                    row_clusters = NULL,
  #                    col_clusters = NULL,
  #                    row_clusters_proba = NULL, #row_clusters_proba_ref,#
  #                    col_clusters_proba = NULL, #col_clusters_proba_ref,#
  #                    X0 = NULL, B0 = NULL,
  #                    min_X0 = 0, max_X0 = 10,
  #                    max_X0B0 = 0.2)



  # res <- one_ZIPLN_simulation(1, "species", simu_params, "Abundance ~ 1 + V1",
  #                             "Abundance ~ 1 + V1", NA)

# Y <- matrix(rep(0, simu_params$n, simu_params$p), nrow = simu_params$n)
# while( (TRUE %in% (rowSums(Y) == 0)) | (TRUE %in% (colSums(Y) == 0)) ){
#   params <- do.call(generate_all_ZIPLN_parameters, simu_params)
#   Y <- simulate_ZIPLN_data(params)
# }
# if(!is.null(params$zi_params)){
#   X <- data.frame(params$X, as.factor(params$zi_params$X0))
#   colnames(X)[[length(colnames(X))]] <- "VZI1"
# }else{X <- params$X}
# simu_data <- prepare_data(Y, X)
# params$Y <- simu_data$Abundance
#
#
# myPLN <- PLNnetwork(as.formula(PLN_formula), simu_data,
#                     control = PLNnetwork_param(penalize_diagonal = FALSE,
#                                                min_ratio = 0.01,
#                                                n_penalties = 20,
#                                                trace = 0))
# PLN_model <- myPLN
#
#
# zi <- ifelse(simu_params$zi_type == "sites", "row",
#              ifelse(simu_params$zi_type == "species", "col", NA) )
# myZIPLN <- ZIPLNnetwork(as.formula(ZIPLN_formula), simu_data, zi = zi,
#                         control = ZIPLNnetwork_param(penalize_diagonal = FALSE,
#                                                      min_ratio = 0.01,
#                                                      n_penalties = 40,
#                                                      trace = 0))

# PLN_formula <- setting$PLN_formula
# ZIPLN_formula <- setting$ZIPLN_formula
# PLN_formula_ZIvar <- setting$PLN_formula_ZIvar
# simu_params = list(n = setting$n,
#                    p = setting$p,
#                    d = 1,
#                    omega_structure = setting$omega_structure,
#                    zi_type = setting$zi_type,
#                    zi_covar_cluster = TRUE,
#                    min_X = 0, max_X = 10, mean_B  = 2, sd_B = 1, XB_max = 70,
#                    n_mode_zi_proba = setting$n_mode_zi_proba,
#                    zi_mode_values = setting$zi_mode_values[[1]],
#                    proba_mode_zi = setting$proba_mode_zi[[1]],
#                    block_values = setting$block_values,
#                    row_clusters = NULL,
#                    col_clusters = NULL,
#                    row_clusters_proba = setting$row_clusters_proba,
#                    col_clusters_proba = setting$col_clusters_proba,
#                    X0 = NULL, B0 = NULL,
#                    min_X0 = 0, max_X0 = 10,
#                    max_X0B0 = 0.2)
#
# res <- one_ZIPLN_simulation(1, zi_config = zi_config, simu_params = simu_params, PLN_formula = PLN_formula,
#                             ZIPLN_formula = ZIPLN_formula,
#                             PLN_formula_ZIvar = PLN_formula_ZIvar)
#
# set.seed(1)

PLN_formula <- "Abundance ~ 0 + V1"
ZIPLN_formula <- "Abundance ~ 0 + V1 | 0 + VZI1" # | 0 + VZI1"
# PLN_formula_ZIvar <- "Abundance ~ 0 " #+ V1 + VZI1

################################ Results plots ################################

# res$zi_strength <- as.numeric(sub(".*_", "", res$zi_config))
# auc_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = AUC, fill = method)) +
#   geom_violin() +
#   facet_grid(zi_strength ~ zi_type) +
#   labs(title = "AUC (zi type ~ zi strength)", x = "", y = "AUC") +
#   theme_minimal() +  # Apply a clean theme
#   scale_fill_viridis_d() +
#   theme(
#     strip.text = element_text(size = 10, face = "bold"),
#     axis.text.x = element_blank(),
#     legend.position = "bottom"
#   )

#################### Network plotting functions ################################
plot_network = function(Omega,
                        type  = "partial_corr",
                        output = "igraph",
                        edge.color      = c("#F8766D", "#00BFC4"),
                        remove.isolated = FALSE,
                        node.labels     = NULL,
                        layout          = layout_in_circle,
                        plot = TRUE) {

  net <- - Omega / tcrossprod(sqrt(diag(Omega))); diag(net) <- 1
  colnames(net) <- unlist(lapply(1:ncol(net), f <- function(i) paste0("species_", i)))
  if (output == "igraph") {

    G <-  graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

    if (!is.null(node.labels)) {
      igraph::V(G)$label <- node.labels
    } else {
      igraph::V(G)$label <- colnames(net)
    }
    ## Nice nodes
    V.deg <- degree(G)/sum(degree(G))
    igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
    igraph::V(G)$size <- V.deg * 100
    igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
    igraph::V(G)$frame.color <- NA
    ## Nice edges
    igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
    if (type == "support"){igraph::E(G)$width <- abs(igraph::E(G)$weight)
    }else{igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)}

    if (remove.isolated) {
      G <- delete.vertices(G, which(degree(G) == 0))
    }
    if (plot) plot(G, layout = layout)
  }
  if (output == "corrplot") {
    if (plot) {
      if (ncol(net) > 100)
        colnames(net) <- rownames(net) <- rep(" ", ncol(net))
      G <- net
      diag(net) <- 0
      corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
    } else  {
      G <- net
    }
  }
  invisible(G)
}

get_partial_corr <- function(Omega){
  p <- nrow(Omega)
  net <- matrix(0, p, p)
  for(i in 1:(p -1)){
    for(j in ((i + 1):p)){
      net[i, j] <- - Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
      net[j, i] <- net[i, j]
    }
  }
  diag(net) <- 1
  return(net)
}

################################## NEGATIVE BINOMIAL DEBUGGING #################
source("negative_binomial_simulate_data.R")
n <- 50 ; p <- 25 ; d <- 1
params <- generate_all_ZIPLN_parameters(n, p, d, add_intercept = TRUE,
                                                    omega_structure = "erdos_renyi",
                                                    zi_type = "species",
                                                    zi_covar_cluster = FALSE,
                                                    min_X = 0, max_X = 1, mean_B  = 2,
                                                    sd_B = 1, XB_max = 70,
                                                    v = 0.3, u = 0.1,
                                                    zi_mode_values = c(0.1),
                                                    proba_mode_zi = c(1),
                                                    block_values = NULL,
                                                    row_clusters = NULL,
                                                    col_clusters = NULL,
                                                    row_clusters_proba = NULL,
                                                    col_clusters_proba = NULL,
                                                    X0 = NULL, B0 = NULL,
                                                    min_X0 = 0, max_X0 = 10,
                                                    max_X0B0 = -0.2)
params$r <- 10
NB_data <- simulate_negative_binomial_data(params)
# pheatmap::pheatmap(NB_data, cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(log(1 + NB_data), cluster_rows = F, cluster_cols = F)
print(median(diag(var(NB_data))))
print(max(diag(var(NB_data))))
print(max(NB_data))

# Poisson_data <- negative_binomial_simulations(params, "poisson")
# # pheatmap::pheatmap(NB_data, cluster_rows = F, cluster_cols = F)
# pheatmap::pheatmap(log(1 + NB_data), cluster_rows = F, cluster_cols = F)
# print(median(diag(var(NB_data))))
# print(max(diag(var(NB_data))))
