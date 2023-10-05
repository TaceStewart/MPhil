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

# Threshold of percent of pre-disturbance cover considered "recovered"
#  Note: for recovery-overlap compounding only
recov_th <- 0.75

# Management benefit (currently % added to mgd reef recovery rate)
mgmt_benefit <- 0.1

# Management constraint (usually 20% of number of reefs in system)
mgmt_constraint <- 0.20

### Sensitivity values ###
sens_recov_yrs <- c(1, 2, 5, 10, 15)
sens_recov_th <- c(0.5, 0.65, 0.75, 0.85, 1)
sens_mgmt_ben <- c(0.02, 0.05, 0.1, 0.15, 0.3)

### Simulations ###
# Run simulations? (Much faster if you don"t)
run_simulations <- TRUE

# Set number of simulations, n_sims
n_sims <- 100

# Set number of sample reefs, num_samples
num_samples <- 1000
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

# Summarise d,r,p for single and compound dist in each sector
sector_df <- data.frame(
  sector = unique(reef_df$management_region),
  num_single = NA,
  ave_prob_s_dist = NA,
  s_dist_ci_lwr = NA,
  s_dist_ci_upr = NA,
  ave_prob_s_impact = NA,
  s_impact_ci_lwr = NA,
  s_impact_ci_upr = NA,
  ave_prob_s_recov = NA,
  s_recov_ci_lwr = NA,
  s_recov_ci_upr = NA,
  num_compound = NA,
  ave_prob_c_dist = NA,
  c_dist_ci_lwr = NA,
  c_dist_ci_upr = NA,
  ave_prob_c_impact = NA,
  c_impact_ci_lwr = NA,
  c_impact_ci_upr = NA,
  ave_prob_c_recov = NA,
  c_recov_ci_lwr = NA,
  c_recov_ci_upr = NA,
  ave_r_s_unmgd = NA,
  ave_r_c_unmgd = NA,
  pr_recov_sing_unmgd = NA,
  pr_recov_comp_unmgd = NA,
  ave_r_s_mgd = NA,
  ave_r_c_mgd = NA,
  pr_recov_sing_mgd = NA,
  pr_recov_comp_mgd = NA
)

# For each section
for (sector in unique(reef_df$management_region)) {
  reef_df_indx <- which(reef_df$management_region == sector)
  sec_df_indx <- which(sector_df$sector == sector)

  # Calculate number of single disturbances
  sector_df$num_single[sec_df_indx] <- reef_df$num_single[reef_df_indx] %>%
    as.numeric() %>%
    sum(na.rm = TRUE)

  # Calculate number of compound disturbances
  sector_df$num_compound[sec_df_indx] <- reef_df$num_comp[reef_df_indx] %>%
    as.numeric() %>%
    sum(na.rm = TRUE)

  # Calculate probability of single disturbance
  sector_df$ave_prob_s_dist[sec_df_indx] <- reef_df$prob_s_dist[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate confidence interval of recovery from single dist (unmanaged)
  ci_s_dist <- reef_df$prob_s_dist[reef_df_indx] %>%
    as.numeric() %>%
    quantile(c(0.025, 0.975), na.rm = TRUE)
  sector_df$s_dist_ci_lwr[sec_df_indx] <- ci_s_dist[1]
  sector_df$s_dist_ci_upr[sec_df_indx] <- ci_s_dist[2]

  # Calculate probability of compound disturbance
  sector_df$ave_prob_c_dist[sec_df_indx] <- reef_df$prob_c_dist[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate confidence interval of recovery from single disturbance
  ci_c_dist <- reef_df$prob_c_dist[reef_df_indx] %>%
    as.numeric() %>%
    quantile(c(0.025, 0.975), na.rm = TRUE)
  sector_df$c_dist_ci_lwr[sec_df_indx] <- ci_c_dist[1]
  sector_df$c_dist_ci_upr[sec_df_indx] <- ci_c_dist[2]

  # Calculate probability of impact from single disturbance
  sector_df$ave_prob_s_impact[sec_df_indx] <- reef_df$prob_s_impact[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate confidence interval of impact from single disturbance
  ci_s_impact <- reef_df$prob_s_impact[reef_df_indx] %>%
    as.numeric() %>%
    quantile(c(0.025, 0.975), na.rm = TRUE)
  sector_df$s_impact_ci_lwr[sec_df_indx] <- ci_s_impact[1]
  sector_df$s_impact_ci_upr[sec_df_indx] <- ci_s_impact[2]

  # Calculate probability of impact from compound disturbance
  sector_df$ave_prob_c_impact[sec_df_indx] <- reef_df$prob_c_impact[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate confidence interval of impact from compound disturbance
  ci_c_impact <- reef_df$prob_c_impact[reef_df_indx] %>%
    as.numeric() %>%
    quantile(c(0.025, 0.975), na.rm = TRUE)
  sector_df$c_impact_ci_lwr[sec_df_indx] <- ci_c_impact[1]
  sector_df$c_impact_ci_upr[sec_df_indx] <- ci_c_impact[2]

  # Calculate probability of recovery from single disturbance (unmanaged)
  sector_df$ave_prob_s_recov[sec_df_indx] <- reef_df$prob_s_recov[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate confidence interval of recovery from single disturbance
  ci_s_recov <- reef_df$prob_s_recov[reef_df_indx] %>%
    as.numeric() %>%
    quantile(c(0.025, 0.975), na.rm = TRUE)
  sector_df$s_recov_ci_lwr[sec_df_indx] <- ci_s_recov[1]
  sector_df$s_recov_ci_upr[sec_df_indx] <- ci_s_recov[2]

  # Calculate probability of recovery from single disturbance (unmanaged)
  sector_df$ave_prob_c_recov[sec_df_indx] <- reef_df$prob_c_recov[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate confidence interval of recovery from single disturbance
  ci_c_recov <- reef_df$prob_c_recov[reef_df_indx] %>%
    as.numeric() %>%
    quantile(c(0.025, 0.975), na.rm = TRUE)
  sector_df$c_recov_ci_lwr[sec_df_indx] <- ci_c_recov[1]
  sector_df$c_recov_ci_upr[sec_df_indx] <- ci_c_recov[2]

  # Calculate probability of recovery from single disturbance (managed)
  ave_r_s_mgd <- reef_df$prob_s_recov[reef_df_indx] %>%
    as.numeric() * (1 + mgmt_benefit) %>%
      min(1)
  sector_df$ave_r_s_mgd[sec_df_indx] <- mean(ave_r_s_mgd, na.rm = TRUE)

  # Calculate probability of recovery from compound disturbance (managed)
  ave_r_c_mgd <- reef_df$prob_c_recov[reef_df_indx] %>%
    as.numeric() * (1 + mgmt_benefit) %>%
      min(1)
  sector_df$ave_r_c_mgd[sec_df_indx] <- mean(ave_r_c_mgd, na.rm = TRUE)

  # Calculate probability of recovery from single disturbance (unmanaged)
  sector_df$ave_r_s_unmgd[sec_df_indx] <- reef_df$prob_s_recov[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate probability of recovery from compound disturbance (unmanaged)
  sector_df$ave_r_c_unmgd[sec_df_indx] <- reef_df$prob_c_recov[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate probability of recovery from single disturbance (unmanaged)
  sector_df$pr_recov_sing_unmgd[sec_df_indx] <- reef_df$pr_recov_sing_unmgd[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate probability of recovery from single disturbance (managed)
  sector_df$pr_recov_sing_mgd[sec_df_indx] <- reef_df$pr_recov_sing_mgd[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate probability of recovery from compound disturbance (unmanaged)
  sector_df$pr_recov_comp_unmgd[sec_df_indx] <- reef_df$pr_recov_comp_unmgd[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)

  # Calculate probability of recovery from compound disturbance (managed)
  sector_df$pr_recov_comp_mgd[sec_df_indx] <- reef_df$pr_recov_comp_mgd[reef_df_indx] %>%
    as.numeric() %>%
    mean(na.rm = TRUE)
}

############################################

########## MGMT AREA OPTIMISATION ##########
# Optimise over management areas
sing_result <- optimiser_single(sector_df, mgmt_constraint) # Single dist sol
comp_result <- optimiser_compound(sector_df, mgmt_constraint) # Compound dist sol

# PRINT RESULTS
reef_s_to_manage_s <- which(get_solution(sing_result, y[i])$value == 1)
prcvd_exp_recov_s <- sum(sector_df$pr_recov_sing_mgd[reef_s_to_manage_s]) +
  sum(sector_df$pr_recov_sing_unmgd[-c(reef_s_to_manage_s)]) # Perceived reefs in recov state
actual_exp_recov_s <- sum(sector_df$pr_recov_comp_mgd[reef_s_to_manage_s]) +
  sum(sector_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_s)]) # Actual reefs in recov state
reef_s_to_manage_c <- which(get_solution(comp_result, y[i])$value == 1)
exp_recov_c <- sum(sector_df$pr_recov_comp_mgd[reef_s_to_manage_c]) +
  sum(sector_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_c)]) # Expected reefs in recov state
print(paste(
  "When only considering single disturbances, it is recommended to manage reef/s",
  reef_s_to_manage_s, "for a perceived", prcvd_exp_recov_s,
  "reefs in a recovered state, but an actual", actual_exp_recov_s,
  "reefs in a recovered state. In the compound disturbance scenario,",
  "it is recommended to manage reef/s", reef_s_to_manage_c,
  "for an expected", exp_recov_c, "reefs in a recovered state."
))
############################################

######## MGMT AREA OPT VISUALISATION #######
# Dumbbell plot
library(ggalt)

sector_df$sector <- as.factor(sector_df$sector)

dumbbell_df <- melt(
  sector_df[, c(
    "sector",
    "pr_recov_comp_unmgd",
    "pr_recov_sing_unmgd"
  )],
  id.vars = "sector"
)
xy <- st_coordinates(st_centroid(sector_boundaries$geometry))
gbr_name_order <- sector_boundaries$AREA_DESCR[order(xy[, "X"], xy[, "Y"],
  decreasing = TRUE
)]

ggplot(
  data = dumbbell_df,
  aes(x = value, y = sector)
) +
  geom_line(color = "azure3") +
  geom_point(aes(color = variable), size = 3) +
  theme(legend.position = "bottom") +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_color_manual(
    name = "",
    labels = c(
      "Compound disturbance included",
      "Singular disturbance only"
    ),
    values = c("steelblue1", "steelblue4")
  ) +
  scale_y_discrete(limits = gbr_name_order[gbr_name_order %in% dumbbell_df$sector]) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(
    x = "Probability the Reef is Recovered",
    y = "GBRMPA Management Area"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = c(.025, .975),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "azure3"),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 5, 4, 3)
  )
if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/Dumbbell_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(), width = 6, height = 5
  )
} else {
  ggsave(
    paste0(
      out_path, "/Dumbbell_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(), width = 6, height = 5
  )
}
############################################

# To Do:
#   make raincloud vis for impact of compounding on prob recovered state
#   Incorporate sample df:
#   Needs calcs for r_single_mgd, r_single_unmgd, r_comp_unmgd, r_comp_mgd,
#                   pr_recov_sing_unmgd, pr_recov_sing_mgd,
#                   pr_recov_comp_unmgd, pr_recov_comp_mgd

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

# PRINT RESULTS
reef_s_to_manage_s <- which(get_solution(sing_result, y[i])$value == 1)
prcvd_exp_recov_s <- sum(sample_reefs_df$pr_recov_sing_mgd[reef_s_to_manage_s]) +
  sum(sample_reefs_df$pr_recov_sing_unmgd[-c(reef_s_to_manage_s)]) # Perceived reefs in recov state
actual_exp_recov_s <- sum(sample_reefs_df$pr_recov_comp_mgd[reef_s_to_manage_s]) +
  sum(sample_reefs_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_s)]) # Actual reefs in recov state
reef_s_to_manage_c <- which(get_solution(comp_result, y[i])$value == 1)
exp_recov_c <- sum(sample_reefs_df$pr_recov_comp_mgd[reef_s_to_manage_c]) +
  sum(sample_reefs_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_c)]) # Expected reefs in recov state
print(paste(
  "When only considering single disturbances, it is recommended to manage reef/s",
  paste(reef_s_to_manage_s, collapse = ", "), "for a perceived",
  prcvd_exp_recov_s, "reefs in a recovered state, but an actual",
  actual_exp_recov_s, "reefs in a recovered state.",
  "In the compound disturbance scenario, it is recommended to manage reef/s",
  paste(reef_s_to_manage_c, collapse = ", "), "for an expected",
  exp_recov_c, "reefs in a recovered state."
))
############################################

############### SIMULATIONS ################

column_names <- c(
  "sector", "point_loc", "point_id", "prob_s_dist", "prob_c_dist", "prob_s_impact",
  "prob_c_impact", "prob_s_recov", "prob_c_recov", "r_single_unmgd", "r_single_mgd",
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
    is_managed_cumul = NA
  ) # column for managed/not (compound)

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

# Visualise the reefs to manage
all_samples <- st_as_sf(all_samples, crs = st_crs(sector_boundaries))
single_vals <- st_transform(all_samples[all_samples$is_managed_single == 1, ],
  crs = st_crs(sector_boundaries)
)

# Create a faceted 2x2 plot of histograms of the probability of single disturbance for each management area
ggplot(all_samples, aes(x = prob_s_dist)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sector, nrow = 2, ncol = 2) +
  labs(x = "Probability of Single Disturbance", y = "Count")

# Create a faceted 2x2 plot of histograms of the probability of compound disturbance for each management area
ggplot(all_samples, aes(x = prob_c_dist)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sector, nrow = 2, ncol = 2) +
  labs(x = "Probability of Compound Disturbance", y = "Count")

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
  geom_hex(single_vals,
    mapping = aes(x_single, y_single),
    binwidth = c(0.5, 0.5)
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
  geom_hex(compound_vals,
    mapping = aes(x_compound, y_compound),
    binwidth = c(0.5, 0.5)
  ) +
  theme(plot.tag = element_text())
compound_plot

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
  theme(plot.tag = element_text()) +
  geom_hex(single_vals,
    mapping = aes(x_single, y_single),
    binwidth = c(0.5, 0.5)
  ) +
  scale_fill_continuous(limits = c(
    min(p$count, q$count),
    max(p$count, q$count)
  )) +
  labs(fill = "Number of \nManaged Reefs")

compound_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    tag = "B"
  ) +
  geom_hex(compound_vals,
    mapping = aes(x_compound, y_compound),
    binwidth = c(0.5, 0.5)
  ) +
  scale_fill_continuous(limits = c(
    min(p$count, q$count),
    max(p$count, q$count)
  )) +
  theme(plot.tag = element_text())

ggarrange(single_plot, compound_plot,
  ncol = 2, nrow = 1,
  common.legend = TRUE,
  legend = "right"
)

if (is_time_based) {
  ggsave(
    paste0(
      out_path, "/",
      n_sims, "Sim_Hexmap_TimeBased",
      recov_yrs, "yr", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(),
    width = 8, height = 5
  )
} else {
  ggsave(
    paste0(
      out_path, "/",
      n_sims, "Sim_Hexmap_RecovBased",
      recov_th, "th", mgmt_benefit, "mgmt.png"
    ),
    plot = last_plot(),
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