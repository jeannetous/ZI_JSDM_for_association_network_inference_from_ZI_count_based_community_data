library(ggplot2)
library(viridis)
library(dplyr)

res_species <- read.csv("ZI_simulations_res/ZI_simulations_res_species_3.csv")
res_sites   <- read.csv("ZI_simulations_res/ZI_simulations_res_sites_3.csv")
res_covar   <- read.csv("ZI_simulations_res/ZI_simulations_res_covar_3.csv")
res1 <- rbind(res_species, res_sites, res_covar)
res <- res1[!grepl("^Error", res1$simu), ]
res <- res[res$n != 25,]
res$zi_strength <- as.numeric(sub(".*_", "", res$zi_config))
res <- res[res$zi_strength <= 4,]
res$neutral <- rep(1, nrow(res))
res.covar <- res[res$zi_type == "covar",]
res[res$zi_type == "covar",]$zi_type <- "covariates"
check <- res %>% group_by(zi_config, n, p, omega_structure, criterion, method) %>% summarise(count = n())

# res <- read.csv("ZIPLN_simus_test_3.csv")
res$AUC <- as.numeric(res$AUC) #; res <- res[!is.na(res$AUC),]
res$f1_score <- as.numeric(res$f1_score)
res$omega_rmse <- as.numeric(res$omega_rmse)
res$precision <- as.numeric(res$precision)
res$fit_rmse <- as.numeric(res$fit_rmse)
res$time <- as.numeric(res$time)
# res_auc <- res[res$criterion == "BIC",]
res <- res[res$method != "ZIPLN_extra_covar",]

res$method_criterion <- paste0(res$method, "_", res$criterion)
res$method_criterion <- sub("_NA$", "", res$method_criterion)
res$method_criterion.2 <- res$method_criterion
res$method_criterion.2 <- sub("_0.9$", "", res$method_criterion.2)
res$method_criterion.2 <- sub("_0.8$", "", res$method_criterion.2)
# res <- res[res$n %in% c(100, 300) & res$method != "neighborhood_selection_network",]
res$method_legend <- res$method
res[res$method_legend == "PLN_ZIvar",]$method_legend <- "PLN with ZI covariate"
res[res$method_legend == "hmsc_ZIvar",]$method_legend <- "Hmsc with ZI covariate"
res[res$method_legend == "gllvm_ZIvar",]$method_legend <- "gllvm with ZI covariate"
res[res$method_legend == "hmsc",]$method_legend <- "Hmsc"

res$method_criterion_legend <- res$method_criterion
res[res$method_criterion_legend == "PLN_StARS_0.8",]$method_criterion_legend <- "PLN (StARS, 0.8)"
res[res$method_criterion_legend == "PLN_StARS_0.9",]$method_criterion_legend <- "PLN (StARS, 0.9)"
res[res$method_criterion_legend == "ZIPLN_StARS_0.8",]$method_criterion_legend <- "ZIPLN (StARS, 0.8)"
res[res$method_criterion_legend == "ZIPLN_StARS_0.9",]$method_criterion_legend <- "ZIPLN (StARS, 0.9)"
res[res$method_criterion_legend == "PLN_BIC",]$method_criterion_legend <- "PLN (BIC)"
res[res$method_criterion_legend == "ZIPLN_BIC",]$method_criterion_legend <- "ZIPLN (BIC)"
res[res$method_criterion_legend == "PLN_ZIvar_StARS_0.8",]$method_criterion_legend <- "PLN with ZI covariate (StARS, 0.8)"
res[res$method_criterion_legend == "PLN_ZIvar_StARS_0.9",]$method_criterion_legend <- "PLN with ZI covariate (StARS, 0.9)"
res[res$method_criterion_legend == "PLN_ZIvar_BIC",]$method_criterion_legend <- "PLN with ZI covariate (BIC)"
res.f1 <- res[!(res$method %in% c("gllvm", "hmsc", "hmsc_ZIvar", "gllvm_ZIvar")),]

res.time <- res[is.na(res$criterion) | res$criterion !="StARS_0.8",]
res.time[!is.na(res.time$method_criterion_legend) & res.time$method_criterion_legend == "PLN (StARS, 0.9)",]$method_criterion_legend <- "PLN (StARS)"
res.time[!is.na(res.time$method_criterion_legend) & res.time$method_criterion_legend == "ZIPLN (StARS, 0.9)",]$method_criterion_legend <- "ZIPLN (StARS)"
res.time[!is.na(res.time$method_criterion_legend) & res.time$method_criterion_legend == "PLN with ZI covariate (StARS, 0.9)",]$method_criterion_legend <- "PLN with ZI covariate (StARS)"
res.time[res.time$method_legend == "Hmsc",]$method_criterion_legend <- "Hmsc"
res.time[res.time$method_legend == "Hmsc with ZI covariate" ,]$method_criterion_legend <- "Hmsc with ZI covariate"
res.time[res.time$method_legend == "gllvm",]$method_criterion_legend <- "gllvm"
res.time[res.time$method_legend == "gllvm with ZI covariate",]$method_criterion_legend <- "gllvm with ZI covariate"

################################ Figures zi_type ~ strength ####################
auc_zi_type_strength <- ggplot(res, aes(x = neutral, y = AUC, fill = method_legend)) +
  geom_boxplot() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "AUC (zero-inflation type ~ zero-inflation strength)",
       x = "", y = "AUC",
       fill = "Method") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 12, face = "bold",),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    text = element_text(family = "TimesNewRomanPSMT",),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12) 
  )
ggsave("ZI_simulations_figures/simulations_auc.png", auc_zi_type_strength, width = 12, height = 10,
       unit = "in", dpi = 800 )


f1score_zi_type_strength <- ggplot(res.f1, aes(x = neutral, y = f1_score, fill = method_criterion_legend)) +
  geom_boxplot() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "F1 score (zero-inflation type ~ zero-inflation strength)", x = "", y = "F1 score",
       fill = "Method") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 12, face = "bold",),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    text = element_text(family = "TimesNewRomanPSMT",),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12) 
  )+
  ylim(0, 1)
ggsave("ZI_simulations_figures/simulations_f1.png", f1score_zi_type_strength, width = 12, height = 10,
       unit = "in", dpi = 800 )



time_zi_type_strength <- ggplot(res.time, aes(x = neutral, y = time, fill = method_criterion_legend)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "Execution time - log scale (zero-inflation type ~ zero-inflation strength)",
       x = "", y = "time (s)", fill = "Method") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 12, face = "bold",),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    text = element_text(family = "TimesNewRomanPSMT",),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12) 
  )+
  scale_y_log10()

ggsave("ZI_simulations_figures/simulations_time.png", time_zi_type_strength, width = 12, height = 10,
       unit = "in", dpi = 800)

omega_rmse_zi_type_strength <- ggplot(res, aes(x = neutral, y = omega_rmse,
                                               fill = method_criterion_legend)) +
  geom_boxplot() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "Omega RMSE, log scale (zero-inflation type ~ zero-inflation strength)",
       x = "", y = "AUC",
       fill = "Method") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 12, face = "bold",),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    text = element_text(family = "TimesNewRomanPSMT",),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12) 
  )
  scale_y_log10()
################################ Figures n ~ p #################################
res.np <- res[res$zi_strength == 3,]
res.np$n <- factor(res.np$n, levels = c(50, 100))
res.np$p <- factor(res.np$p, levels = c(25, 50))

auc_np_species <- ggplot(res.np[res.np$zi_type == "species",], aes(x = neutral, y = AUC, fill = method_legend)) +
  geom_boxplot() +
  facet_grid(n ~ p, labeller = label_both) +
  labs(title = "AUC - species-dependent zero-inflation (n ~ p)",
       x = "", y = "AUC", fill = "Method") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 12, face = "bold",),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    text = element_text(family = "TimesNewRomanPSMT",),
    plot.title = element_text(size = 16),
    legend.position = "none"
  )

auc_np_sites <- ggplot(res.np[res.np$zi_type == "sites",], aes(x = neutral, y = AUC, fill = method_legend)) +
  geom_boxplot() +
  facet_grid(n ~ p, labeller = label_both) +
  labs(title = "AUC - site-dependent zero-inflation (n ~ p)",
       x = "", y = "AUC", fill = "Method") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 12, face = "bold",),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    text = element_text(family = "TimesNewRomanPSMT",),
    plot.title = element_text(size = 16),
    legend.position = "none"
  )

auc_np_covar <- ggplot(res.np[res.np$zi_type == "covariates",], aes(x = neutral, y = AUC, fill = method_legend)) +
  geom_boxplot() +
  facet_grid(n ~ p, labeller = label_both) +
  labs(title = "AUC - covariate-dependent zero-inflation (n ~ p)",
       x = "", y = "AUC", fill = "Method") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 12, face = "bold",),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    text = element_text(family = "TimesNewRomanPSMT",),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12) 
  )

grid.np <- gridExtra::grid.arrange(auc_np_species, auc_np_sites, auc_np_covar)

ggsave("ZI_simulations_figures/simulations_auc_np.png", grid.np, width = 12, height = 15,
       unit = "in", dpi = 800)
################################ DEPRECATED #################################


rmse_omega_zi_type_strength <- ggplot(res, aes(x = neutral, y = omega_rmse, fill = method_criterion)) +
  geom_boxplot() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "Omega RMSE (zi type ~ zi strength)", x = "", y = "Omega RMSE") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    
    panel.grid.major = element_blank()
  ) + 
  scale_y_log10()


rmse_omega_n_p <- ggplot(res, aes(x = n, y = omega_rmse, fill = method)) +
  geom_violin() +
  facet_grid(n ~ p) +
  labs(title = "Omega RMSE (n ~ p)", x = "", y = "Omega RMSE") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) +
  ylim(0, 1)


precision_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = precision, fill = method)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "Precision (zi type ~ zi strength)", x = "", y = "Precision") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )

fit_zi_type_strength <- ggplot(res, aes(x = zi_strength, y = fit_rmse, fill = method_criterion)) +
  geom_violin() +
  facet_grid(zi_strength ~ zi_type) +
  labs(title = "RMSE fit (zi type ~ zi strength)", x = "", y = "RMSE Fit") +
  theme_minimal() +  # Apply a clean theme
  scale_fill_viridis_d() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )
