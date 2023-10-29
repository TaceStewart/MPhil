########### SET UP THE WORKSPACE ###########
# Clear environment
rm(list = ls())

# Clear plots
if (!is.null(dev.list())) dev.off()

# Clear commands
cat("\014")

############################################

############## LOAD LIBARIES ###############
library(sf) # spatial feature handling
library(lwgeom) # spatial feature handling
library(dplyr) # dataframe manipulation
library(varhandle) # handles variables
library(ompr) # Optimization Modeling Package for R
library(ompr.roi) # Optimization Modeling Package for R, ROI solver
library(reshape2) # reshape data
library(nngeo) # k-Nearest Neighbor Join for Spatial Data
library(lubridate) # easy and fast parsing of date-times
library(tidyr) # tidy data
library(ggpubr) # arranging plots
library(ROI.plugin.glpk) # GNU Linear Programming Kit
library(ROI) # R Optimization Infrastructure
library(ggplot2) # creates plots
library(latex2exp) # LaTeX for ggplot
library(patchwork)
library(ggbump) # Delete in not using for opt vis 2
library(GGally)
library(pracma)

# Load disturbance and recovery calculators
source("Stewart_MPhil_single_or_compound.R")
source("Stewart_MPhil_ortiz_r_func.R")
source("Stewart_MPhil_p_calc.R")
source("Stewart_MPhil_optimiser_single.R")
source("Stewart_MPhil_optimiser_compound.R")
source("sampler_v2.R")
source("Stewart_MPhil_dist_finder.R")
source("Stewart_MPhil_analyser.R")

# Set path to QUT Drive
qut_path <- "../OneDrive - Queensland University of Technology"

# Set the mphil path
mphil_path <- paste0(qut_path, "/Documents/MPhil")

# Set the data path
data_path <- paste0(mphil_path, "/Data")

# Set the figure output path
out_path <- paste0(mphil_path, "/Figures")

# Load disturbances and coral obs (.rds made in code/obs_dist_combine.R)
all_reefs_sf <- readRDS(file = paste0(
  data_path,
  "/all_reefs_sf_gaps_filled.rds"
))

# Load the sample reefs dataset
sample_reefs <- readRDS(file = paste0(
  data_path,
  "/sample_reefs.rds"
))

# Load the shapefile
shapefile_path <- paste0(
  data_path,
  "/GBRMPA/Management_Areas_of_the_Great_Barrier",
  "_Reef_Marine_Park/Management_Areas_of_the_Great_",
  "Barrier_Reef_Marine_Park.shp"
)
sector_boundaries <- st_read(shapefile_path,
                             quiet = TRUE
)

############## SET VARIABLES ###############
# Variables
base_recov_th <- 0.75
list_recov_yrs <- seq(1, 15, by = 1)
infer_baseline <- 1 # Should the baseline be inferred if non-existent at start of obs period?
epsilon <- 0.05
baseline_str <- "mean"

# Set GBR MA names in order from north to south
manament_area_names <- c(
  "Far Northern Management Area",
  "Cairns/Cooktown Management Area",
  "Townsville/Whitsunday Management Area",
  "Mackay/Capricorn Management Area"
)

# Set the mphil path
mphil_path <- "../OneDrive - Queensland University of Technology/Documents/MPhil"

# Set the data path
data_path <- paste0(mphil_path, "/Data")

# Set the figure output path
out_path <- paste0(mphil_path, "/Figures/DataChapter/Sensitivity")

# Time based or recovery based compounding?
#  Note: Set to TRUE for time-based compounding or FALSE for recovery-based
is_time_based <- FALSE

# Minimum values considered a "disturbance" in a year
cots_dist <- 1 # cots per tow
cyc_dist <- 10 # hours of 4m wave height
dhw_dist <- 4 # degree heating weeks

# Get the unique reef names
reef_names <- unique(all_reefs_sf$REEF_NAME)

# Initialise counter for unknown recovery times
overall_unknown_count <- 0

# Add columns to reef data frame for storing:
# - whether there are any disturbances each year,
# - single or compound for each disturbance
# - ortiz recovery rate
# - recovery year for each dist
# - recovery time for each dist
all_reefs_sf <- all_reefs_sf %>%
  mutate(
    is_disturbed = (all_reefs_sf$COTS_value >= cots_dist |
                      all_reefs_sf$Hs4MW_value >= cyc_dist |
                      all_reefs_sf$DHW_value >= dhw_dist),
    is_impacted = 0,
    single_or_compound = NA,
    dist_type = NA,
    recov_year = NA,
    recov_time = NA,
    r_given_impact = NA
  )

############################################

########### SENSITIVITY 1 LOOP #############
# comparison of baseline distribution with and without inferring baseline
baseline_infer_bool <- list(0, 1)
reef_df_all <- data.frame(matrix(NA, nrow = 0, ncol = 17))
recov_th <- base_recov_th
baseline_str <- "mean"
i <- 1
for (infer_baseline in baseline_infer_bool) {
  # SINGLE OR COMPOUND DIST
  # Analyse reefs for single or compound disturbances
  # Returns event_counts, reef_df, overall_unknown_count
  sing_comp_list <- analyse_reefs(
    all_reefs_sf, reef_names,
    is_time_based, recov_yrs, recov_th, cots_dist,
    cyc_dist, dhw_dist, infer_baseline, epsilon,
    baseline_str
  )
  
  # Extract the list elements
  event_counts <- sing_comp_list[[1]]
  all_reefs_sf <- sing_comp_list[[2]]
  reef_df <- sing_comp_list[[3]] %>%
    mutate(baseline_infer_bool = infer_baseline)
  
  reef_df_all[(nrow(reef_df_all) + 1):(nrow(reef_df_all) + nrow(reef_df)), 1:ncol(reef_df)] <- reef_df
  i <- i + 1
}
colnames(reef_df_all) <- colnames(reef_df)
############################################

########### DATA CHAPTER SENS 1 ############
# Violin plot of baseline comparison overall and by management area (new)
reef_df_all$baseline_minus_eps <- reef_df_all$baseline_val %>%
  strsplit(", ") %>%
  lapply(first) %>%
  unlist() %>%
  as.numeric()
reef_df_all$initial_bl <- reef_df_all$baseline_minus_eps / (1 - epsilon)

mu <- reef_df_all %>%
  group_by(baseline_infer_bool) %>% 
  summarise(grp.mean = mean(initial_bl, na.rm = TRUE))

# Part A: Overall
dc_sens_1a <- ggplot(data = reef_df_all, 
                     aes(x = initial_bl)) +
  geom_density(aes(y = after_stat(density) * 100,
                   fill = factor(baseline_infer_bool,
                                 levels = c(1, 0))),
               alpha = 0.4,
               color = "black",
               linewidth = 1) +
  geom_vline(aes(xintercept = grp.mean,
                 color = factor(baseline_infer_bool,
                                levels = c(1, 0))), 
             data = mu,
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Initial Baseline Coral Cover (%)", 
       y = "Proportion of Reefs (%)",
       fill = "Baseline Inferred?",
       color = "Baseline Inferred?",
       tag = "A") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_line(linetype = "dashed", "lightgrey")) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"),
                    labels = c("0" = "No", "1" = "Yes")) +
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"),
                     labels = c("0" = "No", "1" = "Yes"))

# Part B: By management area
mu_1b <- reef_df_all %>%
  group_by(baseline_infer_bool, sector) %>% 
  summarise(grp.mean = mean(initial_bl, na.rm = TRUE))
dc_sens_1b <- ggplot(data = reef_df_all, 
                     aes(x = initial_bl)) +
  geom_density(aes(y = after_stat(density) * 100,
                   fill = factor(baseline_infer_bool,
                                 levels = c(1, 0))),
               alpha = 0.4,
               color = "black",
               linewidth = 1) +
  geom_vline(aes(xintercept = grp.mean,
                 color = factor(baseline_infer_bool,
                                levels = c(1, 0))), 
             data = mu_1b,
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Initial Baseline Coral Cover (%)", 
       y = "Proportion of Reefs (%)",
       fill = "Baseline Inferred?",
       color = "Baseline Inferred?",
       tag = "B") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.major.y = element_line(linetype = "dashed", "lightgrey")) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"),
                    labels = c("0" = "No", "1" = "Yes")) +
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"),
                     labels = c("0" = "No", "1" = "Yes")) +
  facet_wrap(~factor(sector, 
                     levels = c("Far Northern Management Area",
                                "Cairns/Cooktown Management Area",
                                "Townsville/Whitsunday Management Area",
                                "Mackay/Capricorn Management Area"),
                     labels = c("Far Northern",
                                "Cairns/Cooktown",
                                "Townsville/Whitsunday",
                                "Mackay/Capricorn")), 
             ncol = 1)

# Combine the plots
dc_sens_1 <- ggarrange(dc_sens_1a, dc_sens_1b, 
                       ncol = 2, common.legend = TRUE, 
                       legend = "bottom")

# Save the plot
ggsave(
  paste0(out_path, "/Sens_1.png"),
  width = 9, height = 6
)
############################################

########### SENSITIVITY 2 LOOP #############
# comparison of baseline distribution with different baseline inference methods
all_reefs_sf <- all_reefs_sf %>%
  mutate(
    is_disturbed = (all_reefs_sf$COTS_value >= cots_dist |
                      all_reefs_sf$Hs4MW_value >= cyc_dist |
                      all_reefs_sf$DHW_value >= dhw_dist),
    is_impacted = 0,
    single_or_compound = NA,
    dist_type = NA,
    recov_year = NA,
    recov_time = NA,
    r_given_impact = NA
  )

baseline_infer_list <- list("min", "mean", "max")
reef_df_all <- data.frame(matrix(NA, nrow = 0, ncol = 17))
recov_th <- base_recov_th
infer_baseline <- 1
i <- 1
for (baseline_str in baseline_infer_list) {
  # SINGLE OR COMPOUND DIST
  # Analyse reefs for single or compound disturbances
  # Returns event_counts, reef_df, overall_unknown_count
  sing_comp_list <- analyse_reefs(
    all_reefs_sf, reef_names,
    is_time_based, recov_yrs, recov_th, cots_dist,
    cyc_dist, dhw_dist, infer_baseline, epsilon,
    baseline_str
  )
  
  # Extract the list elements
  event_counts <- sing_comp_list[[1]]
  all_reefs_sf <- sing_comp_list[[2]]
  reef_df <- sing_comp_list[[3]] %>%
    mutate(baseline_str = baseline_str)
  
  reef_df_all[(nrow(reef_df_all) + 1):(nrow(reef_df_all) + nrow(reef_df)), 1:ncol(reef_df)] <- reef_df
  i <- i + 1
}
colnames(reef_df_all) <- colnames(reef_df)
############################################

########### DATA CHAPTER SENS 2 ############
# Violin plot of baseline comparison overall and by management area (new)
reef_df_all$baseline_minus_eps <- reef_df_all$baseline_val %>%
  strsplit(", ") %>%
  lapply(first) %>%
  unlist() %>%
  as.numeric()
reef_df_all$initial_bl <- reef_df_all$baseline_minus_eps / (1 - epsilon)

mu <- reef_df_all %>%
  group_by(baseline_str) %>% 
  summarise(grp.mean = mean(initial_bl, na.rm = TRUE))

ggplot(data = reef_df_all, 
       aes(x = initial_bl)) +
  geom_density(aes(y = after_stat(density) * 100,
                   fill = as.factor(baseline_str)),
               alpha = 0.4,
               color = "black",
               linewidth = 1) +
  geom_vline(aes(xintercept = grp.mean,
                 color = as.factor(baseline_str)), 
             data = mu,
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Initial Baseline Coral Cover (%)", 
       y = "Proportion of Reefs (%)",
       fill = "Baseline Estimation\nMethod",
       color = "Baseline Estimation\nMethod") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = c(0.9, 0.9),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed", "lightgrey")) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = c("min" = "#F8766D", "mean" = "#00BFC4", "max" = "#619CFF"),
                    labels = c("min" = "Minimum", "mean" = "Mean", "max" = "Maximum")) +
  scale_color_manual(values = c("min" = "#F8766D", "mean" = "#00BFC4", "max" = "#619CFF"),
                     labels = c("min" = "Minimum", "mean" = "Mean", "max" = "Maximum"))

# Save the plot
ggsave(
  paste0(out_path, "/Sens_2.png"),
  width = 8, height = 5
)
############################################

########### SENSITIVITY 3 LOOP #############
baseline_str <- "mean"
list_epsilon <- seq(0, 1, by = 0.05)

reef_df_all <- data.frame(matrix(NA, nrow = 0, ncol = 17))

# Epsilon
i <- 1
for (epsilon in list_epsilon) {
  all_reefs_sf <- all_reefs_sf %>%
    mutate(
      is_disturbed = (all_reefs_sf$COTS_value >= cots_dist |
                        all_reefs_sf$Hs4MW_value >= cyc_dist |
                        all_reefs_sf$DHW_value >= dhw_dist),
      is_impacted = 0,
      single_or_compound = NA,
      dist_type = NA,
      recov_year = NA,
      recov_time = NA,
      r_given_impact = NA
    )
  
  # SINGLE OR COMPOUND DIST
  # Analyse reefs for single or compound disturbances
  # Returns event_counts, reef_df, overall_unknown_count
  sing_comp_list <- analyse_reefs(
    all_reefs_sf, reef_names,
    is_time_based, recov_yrs, recov_th, cots_dist,
    cyc_dist, dhw_dist, infer_baseline, epsilon,
    baseline_str
  )
  
  # Extract the list elements
  event_counts <- sing_comp_list[[1]]
  all_reefs_sf <- sing_comp_list[[2]]
  reef_df <- sing_comp_list[[3]] %>%
    mutate(epsilon_val = epsilon)
  
  reef_df_all[(nrow(reef_df_all) + 1):(nrow(reef_df_all) + nrow(reef_df)), 1:ncol(reef_df)] <- reef_df
  i <- i + 1
}
colnames(reef_df_all) <- colnames(reef_df)
############################################

########### DATA CHAPTER SENS 3 ############
reef_df_all_sens3 <- reef_df_all %>%
  select(c(reef_name, num_total, num_single, num_comp, epsilon_val)) %>%
  melt(id.vars = c("reef_name", "epsilon_val"),
       variable.name = "dist_type",
       value.name = "num_disturbances")
# Linear plot of number of disturbances for each epsilon value (new)
reef_df_all_sens3$epsilon_val <- as.numeric(reef_df_all_sens3$epsilon_val)
reef_df_all_sens3$num_disturbances <- as.numeric(reef_df_all_sens3$num_disturbances)

# get a mean and confidence interval on the number of disturbances for each epsilon value
# and disturbance type
dc_sens_3_df <- reef_df_all_sens3 %>%
  group_by(epsilon_val, dist_type) %>%
  summarise(mean = mean(num_disturbances, na.rm = TRUE),
            sd = sd(num_disturbances, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower = mean - (1.96 * se),
         upper = mean + (1.96 * se),
         sum_sc_means = sum(mean[dist_type != "num_total"], na.rm = TRUE),
         sum_sc_sds = sum(sd[dist_type != "num_total"], na.rm = TRUE),
         sum_sc_n = sum(n[dist_type != "num_total"], na.rm = TRUE),
         sum_sc_se = sum(se[dist_type != "num_total"], na.rm = TRUE),
         sum_sc_lower = sum(lower[dist_type != "num_total"], na.rm = TRUE),
         sum_sc_upper = sum(upper[dist_type != "num_total"], na.rm = TRUE))

dc_sens_3a <- ggplot(data = dc_sens_3_df[dc_sens_3_df$dist_type != "num_total", ]) +
  geom_point(size = 1, 
             aes(x = epsilon_val * 100, 
                 y = mean,
                 color = as.factor(dist_type))) +
  geom_errorbar(aes(x = epsilon_val * 100,
                    ymin = lower,
                    ymax = upper,
                    color = as.factor(dist_type)),
                width = 0.01) +
  geom_point(size = 1, 
             aes(x = epsilon_val * 100, 
                 y = sum_sc_means,
                 color = "num_total_cat")) +
  geom_errorbar(aes(x = epsilon_val * 100,
                    ymin = sum_sc_lower,
                    ymax = sum_sc_upper,
                    color = "num_total_cat"),
                width = 0.01) +
  geom_hline(yintercept = mean(reef_df_all_sens3$num_disturbances[reef_df_all_sens3$dist_type == "num_total"], 
                               na.rm = TRUE),
             color = "grey30",
             linetype = "dotdash") +
  geom_text(aes(x = 27.5,
                y = mean(reef_df_all_sens3$num_disturbances[reef_df_all_sens3$dist_type == "num_total"], 
                         na.rm = TRUE) + 0.025,
                label = "Total Uncategorised Disturbances"),
            color = "grey30",
            size = 10 / .pt,
            hjust = 0.5,
            vjust = 0) +
  labs(x = TeX("Disturbance Threshold ($\\epsilon$) Value (%)"), 
       y = "Average Number of Disturbances",
       color = "Disturbance Type",
       tag = "A") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.major.y = element_line(linetype = "dashed", "lightgrey")) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 3),
                     breaks = seq(0, 3, 0.5)) +
  scale_color_manual(values = c("num_total_cat" = "grey15", "num_single" = "steelblue1", "num_comp" = "steelblue4"),
                     labels = c("num_total_cat" = "Total Categorised", "num_single" = "Single", "num_comp" = "Cumulative"))

# Get the assumed base case (epsilon = 0.05)
base_case <- reef_df_all_sens3[reef_df_all_sens3$epsilon_val == 0.05, ] %>%
  mutate(dist_type = factor(dist_type, levels = c("num_total", "num_single", "num_comp"))) %>%
  filter(dist_type != "num_total")

# Get all others
others <- reef_df_all_sens3[reef_df_all_sens3$epsilon_val != 0.05, ] %>%
  mutate(dist_type = factor(dist_type, levels = c("num_total", "num_single", "num_comp"))) %>%
  filter(dist_type != "num_total")

# Calculate number of reclassified disturbances from base case
others <- others %>%
  group_by(reef_name, dist_type, epsilon_val) %>%
  summarise(num_disturbances = sum(num_disturbances, na.rm = TRUE))
base_case <- base_case %>%
  group_by(reef_name, dist_type) %>%
  summarise(num_disturbances = sum(num_disturbances, na.rm = TRUE))
others <- others %>%
  left_join(base_case, by = c("reef_name", "dist_type")) %>%
  mutate(num_reclassed = num_disturbances.x - num_disturbances.y) %>%
  select(reef_name, dist_type, epsilon_val, num_reclassed)

# Get a mean and confidence interval on the number of reclassified disturbances
# for each epsilon value
dc_sens_3_reclassed <- others %>%
  group_by(epsilon_val, dist_type) %>%
  summarise(mean = mean(num_reclassed, na.rm = TRUE),
            sd = sd(num_reclassed, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower = mean - (1.96 * se),
         upper = mean + (1.96 * se))

# Plot the number of reclassified disturbances from base case on top of dc_sens_3
dc_sens_3b <- ggplot() +
  geom_smooth(data = dc_sens_3_reclassed,
              aes(x = epsilon_val * 100,
                  y = mean,
                  color = as.factor(dist_type)),
              method = "lm",
              se = TRUE) +
  labs(x = TeX("Disturbance Threshold ($\\epsilon$) Value (%)"), 
       y = "Average Number of Reclassified Disturbances",
       color = "Disturbance Type",
       tag = "B") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.major.y = element_line(linetype = "dashed", "lightgrey")) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_color_manual(values = c("num_single" = "steelblue1",
                                "num_comp" = "steelblue4"),
                     labels = c("num_single" = "Single",
                                "num_comp" = "Cumulative"))

# Combine the plots
dc_sens_3 <- ggarrange(dc_sens_3a, dc_sens_3b, 
                       ncol = 2, nrow = 1, 
                       widths = c(1, 1))

# Save the plot
ggsave(
  paste0(out_path, "/Sens_3.png"),
  width = 11, height = 5
)
############################################

########### SENSITIVITY 4 LOOP #############
# Sensitivity to Recovery Threshold
epsilon <- 0.05
list_recov_th <- seq(0, 1, by = 0.05)
reef_df_all <- data.frame(matrix(NA, nrow = 0, ncol = 17))
i <- 1
for (recov_th in list_recov_th) {
  all_reefs_sf <- all_reefs_sf %>%
    mutate(
      is_disturbed = (all_reefs_sf$COTS_value >= cots_dist |
                        all_reefs_sf$Hs4MW_value >= cyc_dist |
                        all_reefs_sf$DHW_value >= dhw_dist),
      is_impacted = 0,
      single_or_compound = NA,
      dist_type = NA,
      recov_year = NA,
      recov_time = NA,
      r_given_impact = NA
    )
  
  # SINGLE OR COMPOUND DIST
  # Analyse reefs for single or compound disturbances
  # Returns event_counts, reef_df, overall_unknown_count
  sing_comp_list <- analyse_reefs(
    all_reefs_sf, reef_names,
    is_time_based, recov_yrs, recov_th, cots_dist,
    cyc_dist, dhw_dist, infer_baseline, epsilon,
    baseline_str
  )
  
  # Extract the list elements
  event_counts <- sing_comp_list[[1]]
  all_reefs_sf <- sing_comp_list[[2]]
  reef_df <- sing_comp_list[[3]] %>%
    mutate(recov_th_val = recov_th)
  
  reef_df_all[(nrow(reef_df_all) + 1):(nrow(reef_df_all) + nrow(reef_df)), 1:ncol(reef_df)] <- reef_df
  i <- i + 1
}
colnames(reef_df_all) <- colnames(reef_df)

############################################

########### DATA CHAPTER SENS 4 ############
# Violin plot of variation in probability of recovery with recovery threshold
reef_df_all$recov_th_val <- as.numeric(reef_df_all$recov_th_val)
reef_df_all$prob_s_recov <- as.numeric(reef_df_all$prob_s_recov)
reef_df_all$prob_c_recov <- as.numeric(reef_df_all$prob_c_recov)

# get a mean and confidence interval on the number of disturbances for each recovery_th value
# and disturbance type
dc_sens_4_df <- reef_df_all %>%
  select(c(reef_name, recov_th_val, prob_s_recov, prob_c_recov)) %>%
  melt(id.vars = c("reef_name", "recov_th_val"),
       variable.name = "dist_type",
       value.name = "prob_recov") %>%
  group_by(recov_th_val, dist_type) %>%
  summarise(mean = mean(as.numeric(prob_recov), na.rm = TRUE),
            sd = sd(as.numeric(prob_recov), na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower = mean - (1.96 * se),
         upper = mean + (1.96 * se))

dc_sens_4 <- ggplot(data = dc_sens_4_df, 
                    aes(x = recov_th_val * 100,
                        y = mean * 100,
                        color = as.factor(dist_type))) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = lower * 100,
                    ymax = upper * 100),
                width = 0.01) +
  geom_smooth(method = "loess",
              se = TRUE) +
  labs(x = TeX("Recovery Threshold ($\\rho$) Value (%)"),
       y = "Probability of Recovery (%)",
       color = "Disturbance Type") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = c(0.05, 0.05),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed", "lightgrey")) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_color_manual(values = c("prob_s_recov" = "steelblue1",
                                "prob_c_recov" = "steelblue4"),
                     labels = c("prob_s_recov" = "Single",
                                "prob_c_recov" = "Cumulative"))

# Save the plot
ggsave(
  paste0(out_path, "/Sens_4.png"),
  width = 6, height = 5
)

############################################