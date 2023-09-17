# Clear environment
rm(list = ls())

# Clear commands
cat("\014")

# Sensitivity analysis
base_recov_yrs <- 2
base_recov_th <- 0.75
base_mgmt_ben <- 0.1
list_recov_yrs <- c(1, 2, 5, 10, 15)
list_recov_th <- c(0.5, 0.65, 0.75, 0.85, 1)
list_mgmt_ben <- c(0.02, 0.05, 0.1, 0.15, 0.2)
manament_area_names <- c("Far Northern Management Area",
                         "Cairns/Cooktown Management Area",
                         "Townsville/Whitsunday Management Area",
                         "Mackay/Capricorn Management Area")

# Recovery-overlap
recov_sens_df <- expand.grid(X=list_recov_th, Y=list_mgmt_ben)
colnames(recov_sens_df) <- c("recov_th", "mgmt_ben")
recov_sens_df <- recov_sens_df %>% mutate(pr_recov = 0,
                                          exp_recov_rfs = 0)
for (par_combo in 1:nrow(recov_sens_df)) {
  recov_th_val <- recov_sens_df$recov_th[par_combo]
  mgmt_ben_val <- recov_sens_df$mgmt_ben[par_combo]
  
  load(paste0("Parameter Run Environments/RecovBased", 
              recov_th_val, "th", 
              mgmt_ben_val, "mgmt.RData"))
  
  recov_sens_df$pr_recov[par_combo] <- mean(all_samples$pr_recov_comp_unmgd, 
                                            na.rm = T)
  
  recov_sens_df$exp_recov_rfs[par_combo] <- exp_recov_c
}

base_pr_recov <- recov_sens_df$pr_recov[recov_sens_df$recov_th == base_recov_th &
                                          recov_sens_df$mgmt_ben == base_mgmt_ben]
pr_recov_df <- filter(recov_sens_df, mgmt_ben == 0.1)
ggplot(pr_recov_df, aes(x = recov_th*100, 
                        y = pr_recov*100)) + 
  geom_line(colour = "azure3", linewidth = 1) +
  geom_point(colour = "steelblue1", size = 3) +
  labs(x = "Recovery Threshold Value (%)",
       y = "Average Probability a Reef is Recovered (%)") +
  theme_classic() +
  scale_x_discrete(limits=list_recov_th*100)

ggsave("../MPhil Thesis/RecovBasedSensitivity_PrRecov.png", 
       plot = last_plot(), 
       width=8, height=5)

base_exp_recov <- recov_sens_df$exp_recov_rfs[recov_sens_df$recov_th == base_recov_th &
                                                recov_sens_df$mgmt_ben == base_mgmt_ben]
ggplot(recov_sens_df, aes(x = as.factor(recov_th*100), 
                          y = as.factor(mgmt_ben*100), 
                          fill = (exp_recov_rfs))) + 
  geom_tile() +
  scale_fill_gradient2(midpoint=base_exp_recov, low="blue", mid="white",
                       high="red", space ="Lab" ) +
  labs(fill = "Average expected \nnumber of reefs",
       x = "Recovery Threshold (%)",
       y = "Management Benefit (%)") +
  theme_classic()

ggsave("../MPhil Thesis/RecovBasedSensitivity_ExpRecov1.png", 
       plot = last_plot(), 
       width=8, height=5)

ggplot(recov_sens_df, aes(x = recov_th*100, 
                          y = exp_recov_rfs)) + 
  geom_line(aes(linetype = as.factor(mgmt_ben*100)), 
            colour = "azure4",
            linewidth = 1)  +
  scale_color_brewer() +
  geom_point(colour = "steelblue1", size = 2.5) +
  labs(linetype = "Management Benefit (%)",
       x = "Recovery Threshold (%)",
       y = "Average Expected Number of Recovered Reefs") +
  theme_classic() +
  scale_x_discrete(limits=list_recov_th*100)

ggsave("../MPhil Thesis/RecovBasedSensitivity_ExpRecov2.png", 
       plot = last_plot(), 
       width=8, height=5)

# Time-overlap
time_sens_df <- expand.grid(X=list_recov_yrs, Y=list_mgmt_ben)
colnames(time_sens_df) <- c("recov_yr", "mgmt_ben")
time_sens_df <- time_sens_df %>% mutate(pr_recov = 0,
                                        exp_recov_rfs = 0)
for (par_combo in 1:nrow(time_sens_df)) {
  recov_yr_val <- time_sens_df$recov_yr[par_combo]
  mgmt_ben_val <- time_sens_df$mgmt_ben[par_combo]
  
  load(paste0("Parameter Run Environments/TimeBased", 
              recov_yr_val, "yrs", 
              mgmt_ben_val, "mgmt.RData"))
  
  time_sens_df$pr_recov[par_combo] <- mean(all_samples$pr_recov_comp_unmgd, 
                                           na.rm = T)
  
  time_sens_df$exp_recov_rfs[par_combo] <- exp_recov_c
}

pr_recov_df <- filter(time_sens_df, mgmt_ben == 0.1)
ggplot(pr_recov_df, aes(x = recov_yr, 
                        y = pr_recov*100)) + 
  geom_line(colour = "azure3", linewidth = 1) +
  geom_point(colour = "steelblue1", size = 3) +
  labs(x = "Recovery Years",
       y = "Average Probability a Reef is Recovered (%)") +
  theme_classic()  + 
  scale_x_discrete(limits=c(1,2,5,10,15))

ggsave("../MPhil Thesis/TimeBasedSensitivity_PrRecov.png", 
       plot = last_plot(), 
       width=8, height=5)

base_exp_recov <- time_sens_df$exp_recov_rfs[time_sens_df$recov_yr == base_recov_yrs &
                                               time_sens_df$mgmt_ben == base_mgmt_ben]
ggplot(time_sens_df, aes(x = as.factor(recov_yr), 
                         y = as.factor(mgmt_ben*100), 
                         fill = (exp_recov_rfs))) + 
  geom_tile() +
  scale_fill_gradient2(midpoint=base_exp_recov, low="blue", mid="white",
                       high="red", space ="Lab" ) +
  labs(fill = "Average Expected \nNumber of Reefs",
       x = "Recovery Years",
       y = "Management Benefit (%)") +
  theme_classic()

ggsave("../MPhil Thesis/TimeBasedSensitivity_ExpRecov1.png", 
       plot = last_plot(), 
       width=8, height=5)

ggplot(time_sens_df, aes(x = recov_yr, 
                         y = exp_recov_rfs)) + 
  geom_line(aes(linetype = as.factor(mgmt_ben*100)), 
            colour = "azure4",
            linewidth = 1)  +
  scale_color_brewer() +
  geom_point(colour = "steelblue1", size = 2.5) +
  labs(linetype = "Management Benefit (%)",
       x = "Recovery Years",
       y = "Average Expected Number of Reefs") +
  theme_classic()

ggsave("../MPhil Thesis/TimeBasedSensitivity_ExpRecov2.png", 
       plot = last_plot(), 
       width=8, height=5)
