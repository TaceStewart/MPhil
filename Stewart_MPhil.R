########### SET UP THE WORKSPACE ###########
# Clear environment
rm(list = ls())

# Clear plots
if (!is.null(dev.list())) dev.off()

# Clear commands
cat("\014")

############################################

##### LOAD LIBRARIES FUNCTIONS & DATA ######
# Load Libraries
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
############################################

###### SET VARIABLES TO CHANGE ON RUN ######
# Should the baseline be inferred if non-existent at start of obs period?
infer_baseline <- 1
epsilon <- 0.05
baseline_str <- "mean"

# Time based or recovery based compounding?
#  Note: Set to TRUE for time-based compounding or FALSE for recovery-based
is_time_based <- FALSE

### Initial values ###
# Number of years assumed to recover from a disturbance
#  Note: for time-overlap compounding only
recov_yrs <- 5

# Threshold percent of pre-disturbance max cover considered "recovered"
#  Note: for recovery-overlap compounding only
recov_th <- 0.75

# Management benefit 
# (currently % added to mgd reef recovery rate)
# NEW (% of annual prob not recovering)
mgmt_benefit <- 0.5

# Management constraint (base: 20% of number of reefs in system)
mgmt_constraint <- 0.2

### Simulations ###
# Run simulations? (Much faster if you don't)
run_simulations <- TRUE

# Set number of simulations, n_sims
n_sims <- 100 #try 10000

# Set number of sample reefs, num_samples
num_samples <- 100
############################################

####### INITIALISE VARIABLES AND DFS #######
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

######### SINGLE OR COMPOUND DIST ##########
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
reef_df <- sing_comp_list[[3]]
overall_unknown_count <- sing_comp_list[[4]]
############################################

########### HISTORICAL ANALYSIS ############
# All single disturbances
single_dists <- all_reefs_sf %>%
  filter(single_or_compound == "Single")

# All reefs with single disturbances
single_reefs <- reef_df %>%
  filter(num_single > 0)

# All compound disturbances
compound_dists <- all_reefs_sf %>%
  filter(single_or_compound == "Compound")

# Reefs with compound disturbances
compound_reefs <- reef_df %>%
  filter(num_comp > 0)

# Game growth rates following single and compound events
single_dists$recov_time <- as.numeric(single_dists$recov_time)
single_reefs$prob_s_recov <- as.numeric(single_reefs$prob_s_recov)
compound_dists$recov_time <- as.numeric(compound_dists$recov_time)
compound_reefs$prob_c_recov <- as.numeric(compound_reefs$prob_c_recov)
#   Average game recovery time after all single disturbances
game_s_dist_mean <- 1 / mean(single_dists$recov_time, na.rm = TRUE)
#   Average game r (single) across reefs
game_s_reef_mean <- mean(single_reefs$prob_s_recov, na.rm = TRUE)
#   Game confidence interval (95 %)
game_s_reef_ci <- quantile(single_reefs$prob_s_recov,
                           c(0.025, 0.975),
                           na.rm = TRUE
)
#   Average game recovery time after all compound disturbances
game_c_dist_mean <- 1 / mean(compound_dists$recov_time,
                             na.rm = TRUE
)
#   Average game r (compound) across reefs
game_c_reef_mean <- mean(compound_reefs$prob_c_recov,
                         na.rm = TRUE
)
#   Game confidence interval (95 %)
game_c_reef_ci <- quantile(compound_reefs$prob_c_recov,
                           c(0.025, 0.975),
                           na.rm = TRUE
)

# Calculate probability of reef being in a recovered state
reef_df <- p_calculator(reef_df, mgmt_benefit)

############################################

################# SAMPLING #################

sample_reefs_df <- samplerv2(
  reef_df, sector_boundaries, sample_reefs, num_samples,
  is_time_based, recov_th, recov_yrs, cots_dist, cyc_dist,
  dhw_dist
)

############################################

########### SAMPLE OPTIMISATION ############

sample_reefs_df <- p_calculator(sample_reefs_df, mgmt_benefit) # Calculate p for the n reefs
sing_result <- optimiser_single(sample_reefs_df, mgmt_constraint) # Single disturbance solution
comp_result <- optimiser_compound(sample_reefs_df, mgmt_constraint) # Compound disturbance solution

# Get the reefs to manage from solutions
sample_reefs_df$is_managed_single <- get_solution(sing_result, y[i])$value
sample_reefs_df$is_managed_cumul <- get_solution(comp_result, y[i])$value

# Check the optimisation makes sense
sample_reefs_df$difference_single <- sample_reefs_df$pr_recov_sing_mgd -
  sample_reefs_df$pr_recov_sing_unmgd
sample_reefs_df$difference_comp <- sample_reefs_df$pr_recov_comp_mgd -
  sample_reefs_df$pr_recov_comp_unmgd

# Transform to sf
sample_reefs_df <- st_as_sf(sample_reefs_df, crs = st_crs(sector_boundaries))

# Get x and y coordinates
sample_reefs_df <- sample_reefs_df %>%
  mutate(
    x = st_coordinates(point_loc)[, "X"],
    y = st_coordinates(point_loc)[, "Y"]
  )

# Make another column for four types of results
for (i in seq_len(nrow(sample_reefs_df))) {
  if (sample_reefs_df$is_managed_single[i] == 0 && sample_reefs_df$is_managed_cumul[i] == 0) {
    sample_reefs_df$scenario_managed[i] <- "Never selected\nfor management"
  } else if (sample_reefs_df$is_managed_single[i] == 1 && sample_reefs_df$is_managed_cumul[i] == 0) {
    sample_reefs_df$scenario_managed[i] <- "General"
  } else if (sample_reefs_df$is_managed_single[i] == 0 && sample_reefs_df$is_managed_cumul[i] == 1) {
    sample_reefs_df$scenario_managed[i] <- "Single and Cumulative"
  } else if (sample_reefs_df$is_managed_single[i] == 1 && sample_reefs_df$is_managed_cumul[i] == 1) {
    sample_reefs_df$scenario_managed[i] <- "Both"
  }
}

cols <- c("General" = "steelblue1", 
          "Single and Cumulative" = "steelblue4",
          "Both" = "#945cee",
          "Never selected\nfor management" = "grey")

# Make visualisation 1: Map of management regions
# Sample reefs that are never chosen to manage in either scenario are grey circles
# Sample reefs that are chosen to manage in single scenario are steel_blue1 coloured circles
# Sample reefs that are chosen to manage in compound scenario are steel_blue4 coloured circles
# Sample reefs that are chosen to manage in both scenarios are steel_blue2 coloured circles
opt_vis_1_1 <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.1) +
  theme_classic() +
  geom_point(
    data = sample_reefs_df,
    aes(x = x, y = y, color = scenario_managed),
    size = 1,
    alpha = 0.75
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  scale_color_manual(
    name = "Scenario Managed",
    values = cols) +
  theme(
    legend.position = c(.025, .025),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "azure3"),
    legend.justification = c("left", "bottom"),
    legend.box.just = "left",
    legend.margin = margin(5, 10, 5, 5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

# Create df for subsetted plot
opt_vis_1_df <- data.frame(
  dist_cons = c("General", "Single and Cumulative", "Single and Cumulative", "None"),
  mgmt_plan = c("General - Perceived", "General - Actual", "Single and Cumulative", "None"),
  benefit_value = c(
    sum(c(sample_reefs_df$pr_recov_sing_mgd[sample_reefs_df$is_managed_single == 1],
          sample_reefs_df$pr_recov_sing_unmgd[sample_reefs_df$is_managed_single == 0])),
    sum(c(sample_reefs_df$pr_recov_comp_mgd[sample_reefs_df$is_managed_single == 1], 
          sample_reefs_df$pr_recov_comp_unmgd[sample_reefs_df$is_managed_single == 0])),
    sum(c(sample_reefs_df$pr_recov_comp_mgd[sample_reefs_df$is_managed_cumul == 1],
          sample_reefs_df$pr_recov_comp_unmgd[sample_reefs_df$is_managed_cumul == 0])),
    sum(sample_reefs_df$pr_recov_comp_unmgd)
  )
)

cols_1_3 <- c("General - Actual" = "steelblue1", 
              "Single and Cumulative" = "steelblue4",
              "General - Perceived" = "lightblue",
              "None" = "grey")
opt_vis_1_3 <- ggplot() +
  geom_bar(
    data = opt_vis_1_df[c(2:4),],
    aes(
      x = reorder(mgmt_plan, benefit_value),
      y = benefit_value,
      fill = mgmt_plan
    ),
    stat = "identity",
    position = position_identity()
  ) +
  geom_text(data = opt_vis_1_df[c(2:4),], 
            aes(x = mgmt_plan,
                y = 15,
                label = round(benefit_value, 2)), 
            vjust = 1,
            col = "white",
            size = 12 / .pt) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"),
    plot.title = element_text(size = 10,
                              hjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_manual(
    name = "",
    values = cols_1_3[c(1,2,4)]
  ) +
  scale_x_discrete(labels = c("None", "General", "Single and\nCumulative")) +
  labs(
    x = "Management Plan",
    title = "Number of Reefs Considered\nRecovered Due to Management"
  )
opt_vis_1_1 + inset_element(opt_vis_1_3, 0.5, 0.6, 0.925, 0.925)

# Save plot
if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/OptVis1_2_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(), width = 8, height = 5
  )
} else {
  ggsave(
    paste0(
      out_path, "/OptVis1_2_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(), width = 8, height = 5
  )
}

# Alter sample_reefs_df to have column of difference in probability of being in a recovered state
# when managed/not
sample_reefs_df <- sample_reefs_df %>%
  mutate(
    diff_prob_recov_s = pr_recov_sing_mgd - pr_recov_sing_unmgd,
    diff_prob_recov_c = pr_recov_comp_mgd - pr_recov_comp_unmgd
  )

# Create dataframe for visualisation 2:
# Point ID | diff_prob_recov_s | diff_prob_recov_c | position_s | position_c | sector
opt_vis_2_df <- data.frame(
  point_id = sample_reefs_df$point_id,
  diff_prob_recov_s = sample_reefs_df$diff_prob_recov_s,
  diff_prob_recov_c = sample_reefs_df$diff_prob_recov_c,
  position_s = NA,
  position_c = NA,
  sector = sample_reefs_df$sector,
  scenario_managed = sample_reefs_df$scenario_managed
)
opt_vis_2_w_ties_s <- rank(length(opt_vis_2_df$diff_prob_recov_s) -
                           opt_vis_2_df$diff_prob_recov_s)
opt_vis_2_w_ties_c <- rank(length(opt_vis_2_df$diff_prob_recov_c) -
                             opt_vis_2_df$diff_prob_recov_c)
a <- data.frame(opt_vis_2_df$diff_prob_recov_s, 
                opt_vis_2_w_ties_s, 
                opt_vis_2_df$diff_prob_recov_c,
                opt_vis_2_w_ties_c)

opt_vis_2_df$position_s <- rank(length(opt_vis_2_df$diff_prob_recov_s) -
                                  opt_vis_2_df$diff_prob_recov_s,
                                ties.method = "random")
opt_vis_2_df$position_c <- rank(length(opt_vis_2_df$diff_prob_recov_c) - 
                                  opt_vis_2_df$diff_prob_recov_c, 
                                ties.method = "random")
opt_vis_2_df$change <- opt_vis_2_df$position_c - opt_vis_2_df$position_s

# Melt the dataframe so that there is one position column and another column for the scenario
opt_vis_2_df <- melt(
  opt_vis_2_df,
  id.vars = c("point_id", "sector", "scenario_managed", "change"),
  measure.vars = c("position_s", "position_c"),
  variable.name = "scenario",
  value.name = "position"
)
opt_vis_2_df <- opt_vis_2_df %>%
  mutate(
    scenario_num = ifelse(scenario == "position_s", 1, 2)
  )

# Make visualisation 2: Bump Chart of the prioritisation of reefs in both management scenarios
# Reefs are ordered by the difference in probability of being in a recovered state when managed/not
cols_vis_2 <- c("Never selected\nfor management" = "grey",
                "General" = "steelblue1",
                "Single and Cumulative" = "steelblue4",
                "Both" = "purple")
opt_vis_2 <- ggplot(opt_vis_2_df, aes(x = scenario, 
                                      y = as.numeric(position), 
                                      color = scenario_managed,
                                      group = point_id)) +
  geom_bump(linewidth = 1.15,
            alpha = 0.65) +
  geom_point(size = 2,
             alpha = 0.75) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "azure3"),
    legend.box = "horizontal",
    legend.margin = margin(c(1,1,1,1), unit = "mm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14,
                                angle = 90),
    axis.text.x = element_text(size = 12,
                               hjust = c(0, 1)),
    axis.text.y = element_blank(),
    plot.tag = element_text(size = 14)
  ) +
  scale_color_manual(values = cols_vis_2) +
  labs(
    x = "Disturbances Considered in Planning Model",
    y = "Reef Priority",
    color = "Model/s Reef is Selected for Management",
    tag = "A"
  ) +
  scale_x_discrete(labels = c("position_s" = "General",
                              "position_c" = "Single and Cumulative"),
                   expand = c(0.02, 0.02),
                   position = "top") +
  scale_y_continuous(expand = c(0.01, 0.01),
                     trans = "reverse") +
  guides(colour = guide_legend(title.position="top", 
                               title.hjust = 0.5,
                               nrow = 1))
opt_vis_2_2 <- ggplot(data = opt_vis_2_df, 
                      mapping = aes(y = change)) +
  geom_histogram(bins = 50,
                 fill = "steelblue4",
                 colour = "black") +
  scale_y_continuous(breaks = seq(round(min(opt_vis_2_df$change), digits = -2), 
                                  round(max(opt_vis_2_df$change), digits = -2), 
                                  by = 10)) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  geom_hline(yintercept = 0,
             colour = "red",
             linewidth = 1) +
  labs(y = "Change in Priority When Incorporating Cumulative Disturbances",
       x = "Number of Reefs",
       tag = "B") +
  geom_text(x = 14.5,
            y = 5,
            label = "Increase in Priority",
            colour = "grey30",
            size = 4) +
  geom_text(x = 14.5,
            y = -4,
            label = "Decrease in Priority",
            colour = "grey30",
            size = 4) +
  theme_pubclean() + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 14))
ggarrange(opt_vis_2, opt_vis_2_2,
          ncol = 2, nrow = 1,
          common.legend = FALSE
)

# Save plot
if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/OptVis2_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(), width = 12, height = 7
  )
} else {
  ggsave(
    paste0(
      out_path, "/OptVis2_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(), width = 12, height = 7
  )
}

############################################

############### SIMULATIONS ################

column_names <- c(
  "sector", "point_loc", "point_id", "num_dist", "num_s_dist", "num_c_dist",
  "prob_s_dist", "prob_c_dist", "prob_s_impact", "prob_c_impact", 
  "prob_s_recov", "prob_c_recov", "r_single_unmgd", "r_single_mgd",
  "r_comp_unmgd", "r_comp_mgd", "pr_recov_sing_unmgd", "pr_recov_sing_mgd",
  "pr_recov_comp_unmgd", "pr_recov_comp_mgd"
)

# Initialise:
all_samples <- matrix(NA, nrow = num_samples * n_sims, ncol = length(column_names)) %>%
  data.frame() # df for all sample reef info
colnames(all_samples) <- column_names
all_samples <- all_samples %>%
  mutate(
    is_managed_single = NA, # column in df for managed/not (single)
    is_managed_cumul = NA, # column for managed/not (compound)
    scenario_managed = NA,
    sim_num = rep(seq(n_sims), each = num_samples) # column for simulation number
  ) 

# For each simulation,
for (sim in 1:n_sims) {
  # Sample num_samples reefs from reef_df
  sample_reefs_df <- samplerv2(
    reef_df, sector_boundaries, sample_reefs, num_samples,
    is_time_based, recov_th, recov_yrs, cots_dist, cyc_dist,
    dhw_dist
  )
  
  # Calculate p for the n reefs
  sample_reefs_df <- p_calculator(sample_reefs_df, mgmt_benefit)
  
  from_row <- (sim - 1) * num_samples + 1
  to_row <- sim * num_samples
  all_samples[from_row:to_row, 1:length(column_names)] <- sample_reefs_df %>%
    subset(select = column_names)
  
  # Calculate single disturbance solution
  sing_result <- optimiser_single(sample_reefs_df, mgmt_constraint)
  
  # Calculate compound disturbance solution
  comp_result <- optimiser_compound(sample_reefs_df, mgmt_constraint)
  
  # Get the reefs to manage from solution
  sol_s <- get_solution(sing_result, y[i])$value
  sol_c <- get_solution(comp_result, y[i])$value
  
  all_samples[from_row:to_row, "is_managed_single"] <- sol_s
  all_samples[from_row:to_row, "is_managed_cumul"] <- sol_c
}

# Make another column for four types of results
for (i in seq_len(nrow(all_samples))) {
  if (all_samples$is_managed_single[i] == 0 && all_samples$is_managed_cumul[i] == 0) {
    all_samples$scenario_managed[i] <- "Never selected\nfor management"
  } else if (all_samples$is_managed_single[i] == 1 && all_samples$is_managed_cumul[i] == 0) {
    all_samples$scenario_managed[i] <- "General"
  } else if (all_samples$is_managed_single[i] == 0 && all_samples$is_managed_cumul[i] == 1) {
    all_samples$scenario_managed[i] <- "Single and Cumulative"
  } else if (all_samples$is_managed_single[i] == 1 && all_samples$is_managed_cumul[i] == 1) {
    all_samples$scenario_managed[i] <- "Both"
  }
}

opt_result_1 <- all_samples %>%
  group_by(scenario_managed) %>%
  summarise(perc = n() / nrow(all_samples))

# Visualise the reefs to manage
all_samples <- st_as_sf(all_samples, crs = st_crs(sector_boundaries))
single_vals <- st_transform(all_samples[all_samples$is_managed_single == 1, ],
                            crs = st_crs(sector_boundaries)
)

# Create a faceted 2x2 plot of histograms of the probability of single disturbance for each management area
ggplot(all_samples, aes(x = prob_s_dist)) +
  geom_histogram(bins = 20) +
  facet_wrap(~sector, nrow = 2, ncol = 2) +
  labs(x = "Probability of Single Disturbance", y = "Count") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 12000, 500))

# Create a faceted 2x2 plot of histograms of the probability of compound disturbance for each management area
ggplot(all_samples, aes(x = prob_c_dist)) +
  geom_histogram(bins = 20) +
  facet_wrap(~sector, nrow = 2, ncol = 2) +
  labs(x = "Probability of Compound Disturbance", y = "Count") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 12000, 500))

# Make visualisation 3: Density of the number of disturbances, split by managed/not managed
# and by management scenario
opt_vis_3 <- ggplot(all_samples, 
                    aes(x = num_dist, 
                        fill = scenario_managed)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    legend.position = c(0.975, 0.975),
    legend.background = element_blank(),
    legend.justification = c("right", "top"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11)
  ) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Number of Disturbances at Reef",
    y = "Density of Reefs",
    fill = "Scenario Managed"
  )

# Save plot
if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/OptVis3_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_3, width = 5, height = 5
  )
} else {
  ggsave(
    paste0(
      out_path, "/OptVis3_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_3, width = 5, height = 5
  )
}

x_single <- st_coordinates(single_vals$point_loc)[, 1]
y_single <- st_coordinates(single_vals$point_loc)[, 2]

single_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    tag = "A"
  ) +
  stat_summary_hex(
    data = single_vals,
    aes(x = x_single, y = y_single, z = is_managed_single),
    fun = ~sum(.x) / n_sims,
    binwidth = c(0.5, 0.5),
    na.rm = TRUE,
    geom = "hex"
  ) +
  theme(plot.tag = element_text())

compound_vals <- st_transform(all_samples[all_samples$is_managed_cumul == 1, ],
                              crs = st_crs(sector_boundaries)
)

x_compound <- st_coordinates(compound_vals$point_loc)[, 1]
y_compound <- st_coordinates(compound_vals$point_loc)[, 2]

compound_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    tag = "B"
  ) +
  stat_summary_hex(
    data = compound_vals,
    aes(x = x_compound, y = y_compound, z = is_managed_cumul),
    fun = ~sum(.x) / n_sims,
    binwidth = c(0.5, 0.5),
    na.rm = TRUE,
    geom = "hex"
  ) +
  theme(plot.tag = element_text())

p <- ggplot_build(single_plot)$data[[2]]
q <- ggplot_build(compound_plot)$data[[2]]

opt_vis_4 <- data.frame(
  sing_or_cumul = c("General", "Single and Cumulative"),
  variable = c("General", "Single and Cumulative"),
  value = c(
    sum(c(all_samples$pr_recov_comp_mgd[all_samples$is_managed_single == 1],
          all_samples$pr_recov_comp_unmgd[all_samples$is_managed_single == 0])) / 
      n_sims,
    sum(c(all_samples$pr_recov_comp_mgd[all_samples$is_managed_cumul == 1],
          all_samples$pr_recov_comp_unmgd[all_samples$is_managed_cumul == 0])) /
      n_sims
  )
)

single_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    tag = "A"
  ) +
  stat_summary_hex(
    data = single_vals,
    aes(x = x_single, y = y_single, z = is_managed_single),
    fun = ~sum(.x) / n_sims,
    binwidth = c(0.5, 0.5),
    na.rm = TRUE,
    geom = "hex"
  ) +
  scale_fill_gradient2(limits = c(min(p$value, q$value),
                                  max(p$value, q$value)),
                       low = "white",
                       mid = "steelblue1",
                       high = "steelblue4") +
  geom_text(aes(x = 146, y = -24,
                label = TeX(
                  paste("Average $E[R_1] = $",
                        round(opt_vis_4$value[opt_vis_4$sing_or_cumul == "General"], 2)), 
                  output = "character")),
            size = 12/.pt,
            parse = TRUE) +
  labs(fill = "Average Number of \nManaged Reefs") +
  theme(plot.tag = element_text())

compound_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    tag = "B"
  ) +
  stat_summary_hex(
    data = compound_vals,
    aes(x = x_compound, y = y_compound, z = is_managed_cumul),
    fun = ~sum(.x) / n_sims,
    binwidth = c(0.5, 0.5),
    na.rm = TRUE,
    geom = "hex"
  ) +
  theme(plot.tag = element_text()) +
  scale_fill_gradient2(limits = c(min(p$value, q$value),
                                  max(p$value, q$value)),
                       low = "white",
                       mid = "steelblue1",
                       high = "steelblue4") +
  geom_text(aes(x = 146, y = -24,
                label = TeX(
                  paste("Average $E[R_2] = $",
                        round(opt_vis_4$value[opt_vis_4$sing_or_cumul == "Single and Cumulative"], 2)), 
                  output = "character")),
            size = 12/.pt,
            parse = TRUE) +
  labs(fill = "Average Number of \nManaged Reefs") +
  theme(plot.tag = element_text())

opt_vis_4 <- ggarrange(single_plot, compound_plot,
                       ncol = 2, nrow = 1,
                       common.legend = TRUE,
                       legend = "right"
)
opt_vis_4
if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/OptVis4_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_4,
    width = 9, height = 5
  )
} else {
  ggsave(
    paste0(
      out_path, "/OptVis4_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_4,
    width = 9, height = 5
  )
}

############################################

############# SAVE ENVIRONMENT #############
if (is_time_based) {
  save.image(paste0(
    data_path,
    "/Parameter Run Environments/",
    "TimeBased", recov_yrs, "yrs",
    mgmt_benefit, "mgmt.RData"
  ))
} else {
  save.image(paste0(
    data_path,
    "/Parameter Run Environments/",
    "RecovBased", recov_th, "th",
    mgmt_benefit, "mgmt.RData"
  ))
}
############################################