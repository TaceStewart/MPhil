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
infer_baseline <- 0
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

# Management benefit (currently % added to mgd reef recovery rate)
mgmt_benefit <- 0.2

# Management constraint (base: 30% of number of reefs in system)
mgmt_constraint <- 0.20

### Sensitivity values ###
sens_recov_yrs <- c(1, 2, 5, 10, 15)
sens_recov_th <- c(0.5, 0.65, 0.75, 0.85, 1)
sens_mgmt_ben <- c(0.02, 0.05, 0.1, 0.15, 0.3)

### Simulations ###
# Run simulations? (Much faster if you don"t)
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
    sample_reefs_df$scenario_managed[i] <- "Single Only"
  } else if (sample_reefs_df$is_managed_single[i] == 0 && sample_reefs_df$is_managed_cumul[i] == 1) {
    sample_reefs_df$scenario_managed[i] <- "Single and Cumulative"
  } else if (sample_reefs_df$is_managed_single[i] == 1 && sample_reefs_df$is_managed_cumul[i] == 1) {
    sample_reefs_df$scenario_managed[i] <- "Both"
  }
}

cols <- c("Single Only" = "steelblue1", 
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
  sing_or_cumul = c("Single Only", "Single Only", "Single and Cumulative", "No Management"),
  variable = c("Single Only - Perceived", "Single Only - Actual", "Single and Cumulative", "No Management"),
  value = c(
    sum(sample_reefs_df$pr_recov_sing_mgd[sample_reefs_df$is_managed_single == 1]) +
      sum(sample_reefs_df$pr_recov_sing_unmgd[sample_reefs_df$is_managed_single == 0]),
    sum(sample_reefs_df$pr_recov_comp_mgd[sample_reefs_df$is_managed_single == 1]) +
      sum(sample_reefs_df$pr_recov_comp_unmgd[sample_reefs_df$is_managed_single == 0]),
    sum(sample_reefs_df$pr_recov_comp_mgd[sample_reefs_df$is_managed_cumul == 1]) +
      sum(sample_reefs_df$pr_recov_comp_unmgd[sample_reefs_df$is_managed_cumul == 0]),
    sum(sample_reefs_df$pr_recov_comp_unmgd)
  )
)

# Add inset plot to vis 1: simple bar chart
# of the difference in expected number of reefs in a recovered state for single and cumulative scenario
# compared to no management
opt_vis_1_df_2 <- data.frame(
  sing_or_cumul = c("Single Only", "Single and Cumulative"),
  variable = c("Single Only", "Single and Cumulative"),
  value = c(
    sum(sample_reefs_df$pr_recov_comp_mgd[sample_reefs_df$is_managed_single == 1]) +
      sum(sample_reefs_df$pr_recov_comp_unmgd[sample_reefs_df$is_managed_single == 0]) -
      sum(sample_reefs_df$pr_recov_comp_unmgd),
    sum(sample_reefs_df$pr_recov_comp_mgd[sample_reefs_df$is_managed_cumul == 1]) +
      sum(sample_reefs_df$pr_recov_comp_unmgd[sample_reefs_df$is_managed_cumul == 0]) -
      sum(sample_reefs_df$pr_recov_comp_unmgd)
  )
)
cols_1_3 <- c("Single Only" = "steelblue1", 
              "Single and Cumulative" = "steelblue4")
opt_vis_1_3 <- ggplot() +
  geom_bar(
    data = opt_vis_1_df_2,
    aes(
      x = sing_or_cumul,
      y = value,
      fill = variable
    ),
    stat = "identity",
    position = position_identity()
  ) +
  geom_text(data = opt_vis_1_df_2, 
            aes(x = sing_or_cumul,
                y = 0.25,
                label = round(value, 2)), 
            vjust = 1,
            col = "white",
            size = 8) +
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
    # legend.position = "bottom",
    # legend.text = element_text(size = 6),
    # legend.box.margin = margin(0, 0, 0, 0),
    # legend.box.spacing = unit(0, "cm")
  ) +
  scale_fill_manual(
    name = "",
    values = cols_1_3
  ) +
  labs(
    x = "Disturbances Considered",
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
    plot = last_plot(), width = 5, height = 5
  )
} else {
  ggsave(
    paste0(
      out_path, "/OptVis1_2_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(), width = 5, height = 5
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
  sector = sample_reefs_df$sector
)
opt_vis_2_df$position_s <- rank(opt_vis_2_df$diff_prob_recov_s, ties.method = "random")
opt_vis_2_df$position_c <- rank(opt_vis_2_df$diff_prob_recov_c, ties.method = "random")

# Melt the dataframe so that there is one position column and another column for the scenario
opt_vis_2_df <- melt(
  opt_vis_2_df,
  id.vars = c("point_id", "sector"),
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
opt_vis_2 <- ggplot(opt_vis_2_df, aes(x = scenario, 
                                      y = as.numeric(position), 
                                      color = gsub("Management Area", "", sector),
                                      group = point_id)) +
  geom_bump(linewidth = 1,
            alpha = 0.5) +
  geom_point(size = 1) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "azure3"),
    legend.justification = c("left", "centre"),
    legend.box.just = "left",
    legend.margin = margin(5, 10, 5, 5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9)
  ) +
  labs(
    x = "Management Scenario",
    y = "Reef Priority",
    color = "GBRMPA Management Area"
  ) +
  scale_y_reverse() +
  scale_x_discrete(labels = c("position_s" = "Single Only",
                              "position_c" = "Single and\nCumulative"),
                   expand = c(0, 0.3))

# Save plot
if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/OptVis2_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_2, width = 6, height = 6
  )
} else {
  ggsave(
    paste0(
      out_path, "/OptVis2_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_2, width = 6, height = 6
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
    all_samples$scenario_managed[i] <- "Single Only"
  } else if (all_samples$is_managed_single[i] == 0 && all_samples$is_managed_cumul[i] == 1) {
    all_samples$scenario_managed[i] <- "Single and Cumulative"
  } else if (all_samples$is_managed_single[i] == 1 && all_samples$is_managed_cumul[i] == 1) {
    all_samples$scenario_managed[i] <- "Both"
  }
}

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
  scale_fill_continuous(limits = c(
    min(p$value, q$value),
    max(p$value, q$value)
  )) +
  labs(fill = "Average Number of \nManaged Reefs")

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
  scale_fill_continuous(limits = c(
    min(p$value, q$value),
    max(p$value, q$value)
  )) +
  theme(plot.tag = element_text())

opt_vis_4 <- ggarrange(single_plot, compound_plot,
                       ncol = 2, nrow = 1,
                       common.legend = TRUE,
                       legend = "right"
)

if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/OptVis4_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_4,
    width = 8, height = 5
  )
} else {
  ggsave(
    paste0(
      out_path, "/OptVis4_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = opt_vis_4,
    width = 8, height = 5
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