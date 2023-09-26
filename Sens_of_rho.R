########### SET UP THE WORKSPACE ###########
# Clear environment
rm(list = ls())

# Clear plots
if(!is.null(dev.list())) dev.off()

# Clear commands
cat("\014")

############################################

##### LOAD LIBRARIES FUNCTIONS & DATA ######
# Load Libraries
library(sf)               # spatial feature handling
library(lwgeom)           # spatial feature handling
library(dplyr)            # dataframe manipulation
library(varhandle)        # handles variables
library(ompr)             # Optimization Modeling Package for R
library(ompr.roi)         # Optimization Modeling Package for R, ROI solver
library(reshape2)         # reshape data
library(nngeo)            # k-Nearest Neighbor Join for Spatial Data
library(lubridate)        # easy and fast parsing of date-times
library(tidyr)            # tidy data
library(ggpubr)           # arranging plots
library(ROI.plugin.glpk)  # GNU Linear Programming Kit
library(ROI)              # R Optimization Infrastructure
library(ggplot2)          # creates plots
library(latex2exp)        # LaTeX for ggplot

# Load disturbance and recovery calculators
source('Stewart_MPhil_single_or_compound.R')
source('Stewart_MPhil_ortiz_r_func.R')
source('Stewart_MPhil_p_calc.R')
source('Stewart_MPhil_Optimiser_Single.R')
source('Stewart_MPhil_Optimiser_compound.R')
source('Stewart_MPhil_Sampler.R')
source('Stewart_MPhil_Dist_Finder.R')

# Set the mphil path 
mphil_path <- "../OneDrive - Queensland University of Technology/Documents/MPhil"

# Set the data path 
data_path <- paste0(mphil_path, "/Data")

# Set the figure output path 
out_path <- paste0(mphil_path, "/Figures/DataChapterOutputs")

# Load disturbances and coral obs (.rds made in code/obs_dist_combine.R)
all_reefs_sf <- readRDS(file = paste0(data_path,
                                      "/all_reefs_sf_gaps_filled.rds"))

# Load the shapefile
shapefile_path <- paste0(data_path,
                         "/GBRMPA/Management_Areas_of_the_Great_Barrier",
                         "_Reef_Marine_Park/Management_Areas_of_the_Great_",
                         "Barrier_Reef_Marine_Park.shp")
sector_boundaries <- st_read(shapefile_path,
                             quiet = TRUE)
map_data <- st_read(shapefile_path)
############################################

rho_values <- seq(0.5, 1, by = 0.05)

# Initialize empty vectors to store results
pr_recov_single_df <- numeric(length(rho_values))
pr_recov_cumulative_df <- numeric(length(rho_values))
i = 1
# start forloop here:
for (recov_th in rho_values) {
  
  ###### SET VARIABLES TO CHANGE ON RUN ######
  # Should the baseline be inferred if non-existent at start of obs period?
  inferBaseline <- 0
  epsilon <- 0.05
  
  # Time based or recovery based compounding?
  #  Note: Set to TRUE for time-based compounding or FALSE for recovery-based
  isTimeBased <- FALSE
  
  ### Initial values ###
  # Number of years assumed to recover from a disturbance 
  #  Note: for time-overlap compounding only
  recov_yrs <- 5
  
  # Threshold for recovery / percent of pre-disturbance cover considered "recovered"
  #  Note: for recovery-overlap compounding only
  #recov_th <- 0.75
  
  # Management benefit (currently % added to mgd reef recovery rate)
  mgmt_benefit <- 0.1
  
  # Management constraint (usually 20% of number of reefs in system)
  mgmt_constraint <- 0.20
  
  ### Sensitivity values ###
  sens_recov_yrs <- c(1, 2, 5, 10, 15)
  sens_recov_th <- c(0.5, 0.65, 0.75, 0.85, 1)
  sens_mgmt_ben <- c(0.02, 0.05, 0.1, 0.15, 0.3)
  
  ### Simulations ###
  # Run simulations? (Much faster if you don't)
  run_simulations <- TRUE
  
  # Set number of simulations, n_sims
  n_sims <- 100
  
  # Set number of reefs, n
  n <- 1
  ############################################
  
  ####### INITIALISE VARIABLES AND DFS #######
  # Minimum values considered a "disturbance" in a year
  cots_dist <- 1 # cots per tow
  cyc_dist <- 10 # hours of 4m wave height
  dhw_dist <- 4  # degree heating weeks
  
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
    mutate(isDisturbed = (all_reefs_sf$COTS_value >= cots_dist | 
                            all_reefs_sf$Hs4MW_value >= cyc_dist |
                            all_reefs_sf$DHW_value >= dhw_dist),
           isImpacted = 0,
           single_or_compound = NA,
           distType = NA,
           recovYear = NA,
           recovTime = NA,
           r_given_impact = NA)
  
  # Initialise df to save reef info in
  reef_df <- data.frame(reef.name = character(),
                        latitude = double(),
                        longitude = double(),
                        sector = character(), 
                        num_single = integer(),
                        prob_s_dist = double(),
                        prob_s_impact = double(),
                        prob_s_recov = double(),
                        num_comp = integer(),
                        prob_c_dist = double(),
                        prob_c_impact = double(),
                        prob_c_recov = double(),
                        reef_unknown = integer())
  
  # Initialise dataframe of disturbance type and counts
  event_counts <- data.frame(H = 0,
                             Cy = 0,
                             Co = 0,
                             HH = 0,
                             HCy = 0,
                             CyH = 0,
                             CyCy = 0,
                             CyCo = 0,
                             CoCy = 0,
                             CoCo = 0,
                             CoH = 0,
                             HCo = 0,
                             Other = 0) # Idea, an "Other" category using paste0 to add remaining compound events
  ############################################
  
  ################# FOR LOOP #################
  # Run code for each reef
  # reef_name <- reef_names[3]
  for (reef_name in reef_names) {
    
    # Extract reef obs for named reef
    obs_by_reef <- all_reefs_sf[all_reefs_sf$REEF_NAME == reef_name,]  
    
    # Separate single and compound disturbances and get recovery times
    if (isTimeBased) {
      obs_by_reef <- single_or_compound(obs_by_reef = obs_by_reef, 
                                        isTimeBased = isTimeBased,
                                        recov_yrs = recov_yrs,
                                        cots_dist = cots_dist,
                                        cyc_dist = cyc_dist,
                                        dhw_dist = dhw_dist,
                                        inferBaseline = inferBaseline,
                                        epsilon = epsilon)
    } else {
      obs_by_reef <- single_or_compound(obs_by_reef = obs_by_reef, 
                                        isTimeBased = isTimeBased,
                                        recov_th = recov_th,
                                        cots_dist = cots_dist,
                                        cyc_dist = cyc_dist,
                                        dhw_dist = dhw_dist,
                                        inferBaseline = inferBaseline,
                                        epsilon = epsilon)
    }
    
    # Calculate number of years observed inclusive
    yrs_obsvd <- obs_by_reef$YEAR[nrow(obs_by_reef)] - obs_by_reef$YEAR[1] + 1
    
    # Calculate disturbance probabilities given occurrences
    # d = (number of events) / (years observed)
    single_dist <- obs_by_reef %>%
      filter(single_or_compound == "Single")
    num_single <- nrow(single_dist)
    prob_s_dist <- num_single/yrs_obsvd
    # Prob Impacted given disturbance
    # P(A|B) = P(A & B)/P(B)
    prob_s_impact <- ifelse(num_single == 0,
                            NA,
                            sum(single_dist$isImpacted & single_dist$isDisturbed,
                                na.rm = T)/
                              sum(single_dist$isDisturbed,
                                  na.rm = T))
    
    prob_s_recov <- ifelse(any(!is.na(single_dist$r_given_impact)),
                           mean(single_dist$r_given_impact, na.rm = T),
                           NA)
    comp_dist <- obs_by_reef %>%
      filter(single_or_compound == "Compound")
    num_comp <- nrow(comp_dist)
    prob_c_dist <- num_comp/yrs_obsvd
    prob_c_impact <- ifelse(num_comp == 0,
                            NA,
                            sum(comp_dist$isImpacted & comp_dist$isDisturbed,
                                na.rm = T)/
                              sum(comp_dist$isDisturbed,
                                  na.rm = T))
    
    prob_c_recov <- ifelse(any(!is.na(comp_dist$r_given_impact)),
                           mean(comp_dist$r_given_impact, na.rm = T),
                           NA)
    
    # Get count of events for the reef and add to total event_counts df
    event_counts_rf <- c(H = sum(obs_by_reef$distType == "Heat Stress", na.rm = TRUE),
                         Cy = sum(obs_by_reef$distType == "Wind Stress", na.rm = TRUE),
                         Co = sum(obs_by_reef$distType == "CoTS Outbreak", na.rm = TRUE),
                         HH = sum(obs_by_reef$distType == "Heat Stress, Heat Stress", na.rm = TRUE),
                         HCy = sum(obs_by_reef$distType == "Heat Stress, Wind Stress", na.rm = TRUE),
                         CyH = sum(obs_by_reef$distType == "Wind Stress, Heat Stress", na.rm = TRUE),
                         CyCy = sum(obs_by_reef$distType == "Wind Stress, Wind Stress", na.rm = TRUE),
                         CyCo = sum(obs_by_reef$distType == "Wind Stress, CoTS Outbreak", na.rm = TRUE),
                         CoCy = sum(obs_by_reef$distType == "CoTS Outbreak, Wind Stress", na.rm = TRUE),
                         CoCo = sum(obs_by_reef$distType == "CoTS Outbreak, CoTS Outbreak", na.rm = TRUE),
                         CoH = sum(obs_by_reef$distType == "CoTS Outbreak, Heat Stress", na.rm = TRUE),
                         HCo = sum(obs_by_reef$distType == "Heat Stress, CoTS Outbreak", na.rm = TRUE),
                         Other = 0)
    event_counts_rf[["Other"]] <- event_counts_rf[["Other"]] + 
      sum(!is.na(obs_by_reef$distType)) - 
      sum(event_counts_rf)
    event_counts <- event_counts + event_counts_rf
    reef_unknown <- sum(grepl("Unknown", obs_by_reef$recovYear), na.rm = TRUE)
    overall_unknown_count <- overall_unknown_count + reef_unknown
    
    ## Get latitude & longitude
    coordinates <- st_coordinates(obs_by_reef$geometry)
    latitude <- coordinates[1,2]
    longitude <- coordinates[1,1]
    
    ## Add all relevant information to data frame 
    new_row <- c(reef_name,                 # reef or reef & site name
                 latitude,                  # latitude of reef/site
                 longitude,                 # longitude of reef/site
                 obs_by_reef$AREA_DESCR[1], # management region
                 num_single,                # no. of singular dists in obs period
                 prob_s_dist,               # probability of single disturbance per year
                 prob_s_impact,             # probability of impact given single disturbance
                 prob_s_recov,              # recovery rate following single dist impact
                 num_comp,                  # no. of cumulative dists in obs period
                 prob_c_dist,               # probability of compound disturbance per year
                 prob_c_impact,             # probability of impact given compound disturbance
                 prob_c_recov,              # recovery rate following compound dist impact
                 reef_unknown)              # number of unknown recovery times
    reef_df[nrow(reef_df) + 1, ] <- new_row
    
    # Replace all_reefs_sf with new info
    all_reefs_sf[all_reefs_sf$REEF_NAME == reef_name,] <- obs_by_reef
  }
  ############################################
  
  ########### HISTORICAL ANALYSIS ############
  # All single disturbances
  single_dists <- all_reefs_sf %>%
    filter(single_or_compound == "Single" &
             isImpacted == 1)
  
  # All compound disturbances
  compound_dists <- all_reefs_sf %>%
    filter(single_or_compound == "Compound" &
             isImpacted == 1)
  
  ############################################
  
  ##########################################
  
  pr_recov_single_df[i] <- mean(as.numeric(reef_df$prob_s_recov), na.rm = T)
  pr_recov_cumulative_df[i] <- mean(as.numeric(reef_df$prob_c_recov), na.rm = T)
  i <- i + 1
}
landscape_dims <- c(8,4)
if (inferBaseline) {
  inferString <- "CS2"
} else {
  inferString <- "CS1"
}

if (isTimeBased) {
  recovString <- paste0("Timebased", recov_yrs, 
                        "yr", mgmt_benefit, "mgmt")
} else {
  recovString <- paste0(recov_th, "th", 
                        mgmt_benefit, "mgmt")
}

# Create a data frame to store the results
results_df <- data.frame(rho = rho_values,
                         num_single = pr_recov_single_df,
                         num_compound = pr_recov_cumulative_df)

# Create a line plot to visualize the results
dumbbell_df <- melt(results_df, 
                    id.vars = "rho")
xy <- st_coordinates(st_centroid(sector_boundaries$geometry))
gbr_name_order <- sector_boundaries$AREA_DESCR[order(xy[,"X"], xy[,"Y"],
                                                     decreasing = TRUE)]
ggplot(data = dumbbell_df, 
       aes(x = rho, y = value, group = rho)) +
  geom_line(color = "azure3") +
  geom_point(aes(color = variable), size = 3) +
  theme(legend.position = "bottom") +
  theme_light() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_color_manual(name = "Disturbance Type", 
                     labels = c("Cumulative", 
                                "Single"),
                     values = c("steelblue1", "steelblue4")) + 
  scale_y_continuous(limits = c(0,1)) +
  labs(x = TeX("Recovery Threshold $(\\rho)$"), 
       y = "Annual Probability of Recovery Following Disturbance") + 
  theme(legend.position = c(.025, .025),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "azure3"),
        legend.justification = c("left", "bottom"),
        legend.box.just = "left")
figname <- "/Sens_of_rho"
ggsave(paste0(out_path, figname, recovString, 
              inferString, ".png"),
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
