########### SET UP THE WORKSPACE ###########
# Clear environment
rm(list = ls())

# Clear plots
if(!is.null(dev.list())) dev.off()

# Clear commands
cat("\014")

############################################

############## LOAD LIBARIES ###############
library(ggplot2)
library(dplyr)

############################################

############## SET VARIABLES ###############
# Variables
base_recov_yrs <- 2
base_recov_th <- 0.75
base_mgmt_ben <- 0.1
list_recov_yrs <- c(1, 2, 5, 10)
list_recov_th <- c(0.5, 0.65, 0.75, 0.85, 1)
list_mgmt_ben <- c(0.02, 0.05, 0.1, 0.15, 0.3)

# Set GBR MA names in order from north to south
manament_area_names <- c("Far Northern Management Area",
                         "Cairns/Cooktown Management Area",
                         "Townsville/Whitsunday Management Area",
                         "Mackay/Capricorn Management Area")

# Set the mphil path 
mphil_path <- "../OneDrive - Queensland University of Technology/Documents/MPhil"

# Set the data path 
data_path <- paste0(mphil_path, "/Data")

# Set the figure output path 
out_path <- paste0(mphil_path, "/Figures/Sensitivity Outputs")

############################################

########### DISTURBANCE HISTORY ############
# For the years considered "compound" in time-overlap, what percentage are: 
#   the same as the recovery-overlap? (true positives + true negatives)
#   false positives?
#   false negatives?
# Need to have recov-overlap base as ground truth
load(paste0(data_path, "/Parameter Run Environments/RecovBased", 
            base_recov_th, "th", 
            base_mgmt_ben, "mgmt.RData"))
ground_truth_df <- data.frame(REEF_NAME = unique(all_reefs_sf$REEF_NAME),
                              single_years = "",
                              compound_years = "")
for (reef in 1:nrow(ground_truth_df)) {
  reef_name <- ground_truth_df$REEF_NAME[reef]
  reef_obs <- all_reefs_sf[all_reefs_sf$REEF_NAME == reef_name,]
  if (any(reef_obs$single_or_compound == "Single", na.rm = T)) {
    single_indices <- which(reef_obs$single_or_compound == "Single")
    for(index in rev(single_indices)) {
      ground_truth_df$single_years[reef] <- reef_obs$recovTime[index] %>% 
        is.na() %>%
        ifelse(max(reef_obs$YEAR), as.integer(reef_obs$recovYear[index]) - 1) %>%
        seq(from = reef_obs$YEAR[index], to = .) %>%
        paste(collapse = ", ") %>%
        paste(ground_truth_df$single_years[reef], sep = ", ")
    }
  } 
  if (any(reef_obs$single_or_compound == "Compound", na.rm = T)) {
    compound_indices <- which(reef_obs$single_or_compound == "Compound")
    for(index in rev(compound_indices)) {
      ground_truth_df$compound_years[reef] <- reef_obs$recovTime[index] %>% 
        is.na() %>%
        ifelse(max(reef_obs$YEAR), as.integer(reef_obs$recovYear[index]) - 1) %>%
        seq(from = reef_obs$YEAR[index], to = .) %>%
        paste(collapse = ", ") %>%
        paste(ground_truth_df$compound_years[reef], sep = ", ")
    }
    
  } 
}

# Need to store years of single and compound dists at each reef in each parcombo
time_sens_parcombos <- expand.grid(X=list_recov_yrs, Y=list_mgmt_ben)
colnames(time_sens_parcombos) <- c("recov_yr", "mgmt_ben")
time_sens_parcombos <- time_sens_parcombos %>%
  mutate(perc_same_single = NA,
         perc_false_single = NA,
         perc_same_comp = NA,
         perc_false_comp = NA)
for (parcombo_row in 1:nrow(time_sens_parcombos)) {
  recov_yr_val <- time_sens_parcombos$recov_yr[parcombo_row]
  mgmt_ben_val <- time_sens_parcombos$mgmt_ben[parcombo_row]
  
  load(paste0(data_path, "/Parameter Run Environments/TimeBased", 
              recov_yr_val, "yrs", 
              mgmt_ben_val, "mgmt.RData"))
  same_single_vec <- c()
  false_single_vec <- c()
  same_compound_vec <- c()
  false_compound_vec <- c()
  # Create df similar to ground truth
  for (reef in 1:nrow(ground_truth_df)) {
    reef_name <- ground_truth_df$REEF_NAME[reef]
    reef_obs <- all_reefs_sf[all_reefs_sf$REEF_NAME == reef_name,]
    reef_single_yrs <- c()
    if (any(reef_obs$single_or_compound == "Single", na.rm = T)) {
      single_indices <- which(reef_obs$single_or_compound == "Single")
      for(index in rev(single_indices)) {
        reef_single_yrs <- reef_obs$recovTime[index] %>% 
          is.na() %>%
          ifelse(max(reef_obs$YEAR), as.integer(reef_obs$recovYear[index]) - 1) %>%
          seq(from = reef_obs$YEAR[index], to = .) %>%
          append(reef_single_yrs)
      }
    } 
    if (any(reef_obs$single_or_compound == "Compound", na.rm = T)) {
      compound_indices <- which(reef_obs$single_or_compound == "Compound")
      reef_compound_yrs <- c()
      for(index in rev(compound_indices)) {
        reef_compound_yrs <- reef_obs$recovTime[index] %>% 
          is.na() %>%
          ifelse(max(reef_obs$YEAR), as.integer(reef_obs$recovYear[index]) - 1) %>%
          seq(from = reef_obs$YEAR[index], to = .)  %>%
          append(reef_compound_yrs)
      }
    } 
    
    # Evaluate % same single and compound
    ground_truth_reef_single <- ground_truth_df$single_years[reef] %>%
      strsplit(", ") %>%
      unlist() %>%
      as.numeric()
    ground_truth_reef_comp <- ground_truth_df$compound_years[reef] %>%
      strsplit(", ") %>%
      unlist() %>%
      as.numeric()
    
    # Calc % of the parcombo reef single years are in gt single years
    if (length(ground_truth_reef_single) > 0) {
      perc_same_single <- sum(reef_single_yrs %in% ground_truth_reef_single)/
        length(ground_truth_reef_single)
      same_single_vec <- ifelse(is.nan(perc_same_single),
                                0, perc_same_single) %>%
        append(same_single_vec)
    } else {
      same_single_vec <- ifelse(length(reef_single_yrs) > 0,
                                0, 1) %>%
        append(same_single_vec)
    }
    
    # Calc % of the parcombo reef single years are not in gt single years
    if (length(reef_single_yrs) > 0) {
      perc_false_single <- sum(!(reef_single_yrs %in% ground_truth_reef_single))/
        length(reef_single_yrs)
      false_single_vec <- ifelse(is.nan(perc_same_single),
                                 0, perc_same_single) %>%
        append(false_single_vec)
    } else {
      false_single_vec <- ifelse(length(ground_truth_reef_single) > 0,
                                 0,1) %>%
        append(false_single_vec)
    }
    
    # Calc % of the parcombo reef compound years are in gt compound years
    if (length(ground_truth_reef_comp) > 0) {
      perc_same_comp <- sum(reef_compound_yrs %in% ground_truth_reef_comp)/
        length(ground_truth_reef_comp)
      same_compound_vec <- ifelse(is.nan(perc_same_comp),
                                  0, perc_same_comp) %>%
        append(same_compound_vec)
    } else {
      same_compound_vec <- ifelse(length(reef_compound_yrs) > 0,
                                  0,1) %>%
        append(same_compound_vec)
    }
    
    # Calc % of the parcombo reef compound years are not in gt compound years
    if (length(reef_compound_yrs) > 0) {
      perc_false_comp <- sum(!(reef_compound_yrs %in% ground_truth_reef_comp))/
        length(reef_compound_yrs)
      false_compound_vec <- ifelse(is.nan(perc_false_comp),
                                   0, perc_false_comp) %>%
        append(false_compound_vec)
    } else {
      false_compound_vec <- ifelse(length(ground_truth_reef_comp) > 0,
                                   0,1) %>%
        append(false_compound_vec)
    }
  }
  
  time_sens_parcombos$perc_same_single[parcombo_row] <- mean(same_single_vec)
  time_sens_parcombos$perc_false_single[parcombo_row] <- mean(false_single_vec)
  time_sens_parcombos$perc_same_comp[parcombo_row] <- mean(same_compound_vec)
  time_sens_parcombos$perc_false_comp[parcombo_row] <- mean(false_compound_vec)
}
time_sens_plot_df <- melt(time_sens_parcombos[time_sens_parcombos$mgmt_ben == 0.1,], 
                id = c("recov_yr", "mgmt_ben"))
ggplot(time_sens_plot_df,
       aes(x = variable, 
           y = value, 
           fill = as.factor(recov_yr))) +
  geom_bar(stat = 'identity',
           position = 'dodge') +
  theme_classic() + 
  scale_fill_manual(values=brewer.pal(9, "Blues")[3:6]) +
  labs(x = "Metric Measured Against Recovery-Overlap",
       y = "Percent",
       fill  = "Estimated\nRecovery\nTime (Years)") +
  scale_x_discrete(labels = c("Same Single\nDisturbances",
                              "False Single\nDisturbances",
                              "Same Compound\nDisturbances",
                              "False Compound\nDisturbances"))

ggsave(paste0(out_path, "/DHSensitivity_RecoveryTime.png"), 
       plot = last_plot(), 
       width=8, height=5)


############################################

############# OPTIMAL PLANNING ##############
# Recovery-overlap
recov_sens_df <- expand.grid(X=list_recov_th, Y=list_mgmt_ben)
colnames(recov_sens_df) <- c("recov_th", "mgmt_ben")
recov_sens_df <- recov_sens_df %>% mutate(pr_recov = 0,
                                          exp_recov_rfs = 0)
for (par_combo in 1:nrow(recov_sens_df)) {
  recov_th_val <- recov_sens_df$recov_th[par_combo]
  mgmt_ben_val <- recov_sens_df$mgmt_ben[par_combo]
  
  load(paste0(data_path, "/Parameter Run Environments/RecovBased", 
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

ggsave(paste0(out_path, "/RecovBasedSensitivity_PrRecov.png"), 
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

ggsave(paste0(out_path, "/RecovBasedSensitivity_ExpRecov1.png"), 
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

ggsave(paste0(out_path, "/RecovBasedSensitivity_ExpRecov2.png"), 
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
  
  load(paste0(data_path, "/Parameter Run Environments/TimeBased", 
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

ggsave(paste0(data_path, "/TimeBasedSensitivity_PrRecov.png"), 
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

ggsave(paste0(data_path, "/TimeBasedSensitivity_ExpRecov1.png"), 
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

ggsave(paste0(data_path, "/TimeBasedSensitivity_ExpRecov2.png"), 
       plot = last_plot(), 
       width=8, height=5)
############################################



