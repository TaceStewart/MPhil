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
source('Stewart_MPhil_Optimiser_Compound.R')
source('Stewart_MPhil_Sampler.R')
source('Stewart_MPhil_Dist_Finder.R')

# Set the mphil path 
mphil_path <- "../OneDrive - Queensland University of Technology/Documents/MPhil"

# Set the data path 
data_path <- paste0(mphil_path, "/Data")

# Set the figure output path 
out_path <- paste0(mphil_path, "/Figures/Data Chapter Outputs")

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
############################################

###### SET VARIABLES TO CHANGE ON RUN ######
# Time based or recovery based compounding?
#  Note: Set to TRUE for time-based compounding or FALSE for recovery-based
isTimeBased <- FALSE

### Initial values ###
# Number of years assumed to recover from a disturbance 
#  Note: for time-overlap compounding only
recov_yrs <- 5

# Threshold for recovery / percent of pre-disturbance cover considered "recovered"
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
# Run simulations? (Much faster if you don't)
run_simulations <- TRUE

# Set number of simulations, n_sims
n_sims <- 100

# Set number of reefs, n
n <- 100
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
         single_or_compound = NA,
         distType = NA,
         ortizRecov = NA,
         recovYear = NA,
         recovTime = NA)

# Initialise df to save reef info in
reef_df <- data.frame(reef.name = character(),
                      num_single = integer(),
                      d_single = double(),
                      num_comp = integer(),
                      d_comp = double(),
                      r_ortiz_s = double(),
                      r_ortiz_c = double(),
                      r_single = double(),
                      r_comp = double(),
                      latitude = double(),
                      longitude = double(),
                      sector = character(), 
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
# reef_name <- reef_names[1]
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
                                      dhw_dist = dhw_dist)
  } else {
    obs_by_reef <- single_or_compound(obs_by_reef = obs_by_reef, 
                                      isTimeBased = isTimeBased,
                                      recov_th = recov_th,
                                      cots_dist = cots_dist,
                                      cyc_dist = cyc_dist,
                                      dhw_dist = dhw_dist)
  }
  
  # Calculate Ortiz r for single and compound dists for reef
  obs_by_reef <- ortiz_r_func(obs_by_reef, cots_dist = cots_dist,
                              cyc_dist = cyc_dist, dhw_dist = dhw_dist)
  
  # Some recovery rates will be < 0 because of the threshold. 
  # In these cases, let it be a low r instead
  obs_by_reef$ortizRecov[obs_by_reef$ortizRecov < 0] <- 0.001
  
  # Calculate number of years observed inclusive
  yrs_obsvd <- obs_by_reef$YEAR[nrow(obs_by_reef)] - obs_by_reef$YEAR[1] + 1
  
  # Calculate disturbance probabilities given occurrences
  # d = (number of events) / (years observed)
  d_single <- ifelse(sum(obs_by_reef$single_or_compound == "Single", na.rm = TRUE) == 0,
                     0,
                     sum(obs_by_reef$single_or_compound == "Single", na.rm = TRUE)/yrs_obsvd)
  d_comp <- ifelse(sum(obs_by_reef$single_or_compound == "Compound", na.rm = TRUE) == 0,
                   0,
                   sum(obs_by_reef$single_or_compound == "Compound", na.rm = TRUE)/yrs_obsvd)
  
  # Calculate average r for single and compound events, 
  # r = 1/(ave years to recover)
  single_dist <- obs_by_reef[obs_by_reef$single_or_compound == "Single",]
  single_dist <- single_dist[!is.na(single_dist$recovYear),]
  num_single <- nrow(single_dist)
  if (num_single > 0 & 
      any(!is.na(single_dist$recovTime)) & 
      (any(!grepl("Unknown", single_dist$ortizRecov)) |
       any(!grepl("Unknown", single_dist$recovTime)))) {
    single_dist <- single_dist[!grepl("Unknown", single_dist$ortizRecov),]
    r_ortiz_s <- mean(as.numeric(single_dist$ortizRecov[check.numeric(single_dist$ortizRecov, na.rm = TRUE) &
                                                          single_dist$single_or_compound == "Single"]), 
                      na.rm = TRUE)
    r_single <- 1/mean(as.numeric(single_dist$recovTime[check.numeric(single_dist$recovTime, na.rm = TRUE) &
                                                          single_dist$single_or_compound == "Single"]), 
                       na.rm = TRUE)
  } else {
    r_ortiz_s <- NA
    r_single <- NA
  }
  comp_dist <- obs_by_reef[obs_by_reef$single_or_compound == "Compound",]
  comp_dist <- comp_dist[!is.na(comp_dist$recovYear),]
  num_comp <- nrow(comp_dist)
  if (num_single > 0 & 
      any(!is.na(comp_dist$recovTime)) & 
      (any(!grepl("Unknown", comp_dist$ortizRecov)) |
       any(!grepl("Unknown", comp_dist$recovTime)))) {
    comp_dist <- comp_dist[!grepl("Unknown", comp_dist$ortizRecov),]
    r_ortiz_c <- mean(as.numeric(comp_dist$ortizRecov[check.numeric(comp_dist$ortizRecov, na.rm = TRUE) &
                                                        comp_dist$single_or_compound == "Compound"]), 
                      na.rm = TRUE)
    r_comp <- 1/mean(as.numeric(comp_dist$recovTime[check.numeric(comp_dist$recovTime, na.rm = TRUE) &
                                                      comp_dist$single_or_compound == "Compound"]), 
                     na.rm = TRUE)
  } else {
    num_comp <- 0
    r_ortiz_c <- NA
    r_comp <- NA
  }
  
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
  latitude <- coordinates[2]
  longitude <- coordinates[1]
  
  ## Add all relevant information to data frame 
  new_row <- c(reef_name,         # reef or reef & site name
               num_single,        # no. of singular dists in obs period
               d_single,          # probability of single cyclones per year
               num_comp,          # no. of cumulative dists in obs period
               d_comp,            # probability of compound cyclones per year
               r_ortiz_s,         # ortiz recovery rate (single)
               r_ortiz_c,         # ortiz recovery rate (compound)
               r_single,          # recovery rate following single dist
               r_comp,            # recovery rate following compound dist
               latitude,          # latitude of reef/site
               longitude,         # longitude of reef/site
               obs_by_reef$AREA_DESCR[1], # management region
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
           !grepl("Unknown", recovYear))

all_reefs_sf[all_reefs_sf$single_or_compound == "Single",]

# Reefs with single disturbances
single_reefs <- reef_df[!is.na(reef_df$r_single),] 

# All compound disturbances
compound_dists <- all_reefs_sf[all_reefs_sf$single_or_compound == "Compound" & 
                                 check.numeric(all_reefs_sf$ortizRecov),]
compound_dists <- compound_dists[!is.na(compound_dists$recovYear),] 

# Reefs with compound disturbances
compound_reefs <- reef_df[!is.na(reef_df$r_comp),] 

# Ortiz growth rates following single and compound events
single_dists$ortizRecov <- as.numeric(single_dists$ortizRecov)
single_reefs$r_ortiz_s <- as.numeric(single_reefs$r_ortiz_s)
compound_dists$ortizRecov <- as.numeric(compound_dists$ortizRecov)
compound_reefs$r_ortiz_c <- as.numeric(compound_reefs$r_ortiz_c)
#   Average ortiz r after all single disturbances 
ortiz_s_dist_mean <- mean(single_dists$ortizRecov, na.rm = TRUE) 
#   Average ortiz r (single) across reefs 
ortiz_s_reef_mean <- mean(single_reefs$r_ortiz_s, na.rm = TRUE) 
#   Average ortiz r after all compound disturbances 
ortiz_c_dist_mean <- mean(compound_dists$ortizRecov, na.rm = TRUE) 
#   Average ortiz r (compound) across reefs 
ortiz_c_reef_mean <- mean(compound_reefs$r_ortiz_c, na.rm = TRUE) 

# Game growth rates following single and compound events
single_dists$recovTime <- as.numeric(single_dists$recovTime)
single_reefs$r_single <- as.numeric(single_reefs$r_single)
compound_dists$recovTime <- as.numeric(compound_dists$recovTime)
compound_reefs$r_comp <- as.numeric(compound_reefs$r_comp)
#   Average game recovery time after all single disturbances 
game_s_dist_mean <- 1/mean(single_dists$recovTime, na.rm = TRUE)
#   Average game r (single) across reefs 
game_s_reef_mean <- mean(single_reefs$r_single, na.rm = TRUE) 
#   Game confidence interval (95 %)
game_s_reef_CI <- quantile(single_reefs$r_single, c(0.025, 0.975), na.rm = TRUE)
#   Average game recovery time after all compound disturbances 
game_c_dist_mean <- 1/mean(compound_dists$recovTime, na.rm = TRUE) 
#   Average game r (compound) across reefs 
game_c_reef_mean <- mean(compound_reefs$r_comp, na.rm = TRUE)
#   Game confidence interval (95 %)
game_c_reef_CI <- quantile(compound_reefs$r_comp, c(0.025, 0.975), na.rm = TRUE)

# Data frame to compare recovery following all disturbance events
eventNames <- c("Wind Stress", "Heat Stress", "CoTS Outbreak",
                "Wind Stress, Wind Stress", "Wind Stress, Heat Stress",
                "Wind Stress, CoTS Outbreak", "Heat Stress, Wind Stress",
                "Heat Stress, Heat Stress", "Heat Stress, CoTS Outbreak",
                "CoTS Outbreak, Wind Stress", "CoTS Outbreak, Heat Stress", 
                "CoTS Outbreak, CoTS Outbreak")
recov_df <- data.frame(eventName = eventNames,
                       distType = "",
                       numDist = 0,
                       ortiz_r = 0,
                       r = 0)
for (row in 1:nrow(recov_df)) {
  recov_df$distType[row] <- ifelse(grepl(",", recov_df$eventName[row]), "Compound", "Single")
  recov_df$numDist[row] <- ifelse(recov_df$distType[row] == "Single", 
                                  sum(single_dists$distType == recov_df$eventName[row]),
                                  sum(compound_dists$distType == recov_df$eventName[row]))
  recov_df$ortiz_r[row] <- ifelse(recov_df$distType[row] == "Single",
                                  mean(single_dists$ortizRecov[
                                    single_dists$distType == recov_df$eventName[row]], 
                                    na.rm = T),
                                  mean(compound_dists$ortizRecov[
                                    compound_dists$distType == recov_df$eventName[row]], 
                                    na.rm = T))
  recov_df$r[row] <- ifelse(recov_df$distType[row] == "Single",
                            1/mean(single_dists$recovTime[
                              single_dists$distType == recov_df$eventName[row]], 
                              na.rm = T),
                            1/mean(compound_dists$recovTime[
                              compound_dists$distType == recov_df$eventName[row]], 
                              na.rm = T))
}

# Calculate probability of reef being in a recovered state
reef_df <- p_calculator(reef_df, mgmt_benefit)

# Summarise d,r,p for single and compound dist in each sector
sector_df <- data.frame(sector = unique(reef_df$sector),
                        num_single = NA,
                        num_compound = NA,
                        ave_d_s = NA,
                        ave_d_c = NA,
                        ave_r_s_unmgd = NA,
                        ci_r_s_unmgd_lwr = NA,
                        ci_r_s_unmgd_upr = NA,
                        ave_r_s_mgd = NA,
                        ci_r_s_mgd_lwr = NA,
                        ci_r_s_mgd_upr = NA,
                        ave_r_c_unmgd = NA,
                        ci_r_c_unmgd_lwr = NA,
                        ci_r_c_unmgd_upr = NA,
                        ave_r_c_mgd = NA,
                        ci_r_c_mgd_lwr = NA,
                        ci_r_c_mgd_upr = NA,
                        pr_recov_sing_mgd = NA,
                        pr_recov_sing_unmgd = NA,
                        pr_recov_comp_mgd = NA,
                        pr_recov_comp_unmgd = NA)
# For each section,
# sector <- unique(reef_df$sector)[1]
for (sector in unique(reef_df$sector)) {
  reef_df_indx <- which(reef_df$sector == sector)
  sec_df_indx <- which(sector_df$sector == sector)
  
  # Calculate number of single disturbances
  sector_df$num_single[sec_df_indx] <- sum(as.numeric(reef_df$num_single[reef_df_indx]),
                                           na.rm=TRUE)
  
  # Calculate number of compound disturbances
  sector_df$num_compound[sec_df_indx] <- sum(as.numeric(reef_df$num_comp[reef_df_indx]),
                                             na.rm=TRUE)
  
  # Calculate probability of single disturbance
  sector_df$ave_d_s[sec_df_indx] <- mean(as.numeric(reef_df$d_single[reef_df_indx]),
                                         na.rm=TRUE)
  
  # Calculate probability of compound disturbance
  sector_df$ave_d_c[sec_df_indx] <- mean(as.numeric(reef_df$d_comp[reef_df_indx]),
                                         na.rm=TRUE)
  
  # Calculate probability of recovery from single disturbance (unmanaged)
  sector_df$ave_r_s_unmgd[sec_df_indx] <- mean(as.numeric(reef_df$r_single_unmgd[reef_df_indx]),
                                               na.rm=TRUE)
  
  # Calculate confidence interval of recovery from single disturbance (unmanaged)
  CI <- quantile(as.numeric(reef_df$r_single_unmgd[reef_df_indx]), c(0.025, 0.975), na.rm = TRUE)
  sector_df$ci_r_s_unmgd_lwr[sec_df_indx] <- CI[1]
  sector_df$ci_r_s_unmgd_upr[sec_df_indx] <- CI[2]
  
  # Calculate probability of recovery from single disturbance (managed)
  ave_r_s_mgd <- reef_df$r_single_mgd[reef_df_indx] %>%
    as.numeric() * (1 + mgmt_benefit) %>%
    min(1) 
  sector_df$ave_r_s_mgd[sec_df_indx] <- mean(ave_r_s_mgd, na.rm=TRUE) 
  
  # Calculate confidence interval of recovery from single disturbance (managed)
  CI <- quantile(ave_r_s_mgd, c(0.025, 0.975), na.rm = TRUE)
  sector_df$ci_r_s_mgd_lwr[sec_df_indx] <- CI[1]
  sector_df$ci_r_s_mgd_upr[sec_df_indx] <- CI[2] 
  
  # Calculate probability of recovery from compound disturbance (unmanaged)
  sector_df$ave_r_c_unmgd[sec_df_indx] <- mean(as.numeric(reef_df$r_comp[reef_df_indx]),
                                               na.rm=TRUE)
  
  # Calculate confidence interval of recovery from compound disturbance (unmanaged)
  CI <- quantile(as.numeric(reef_df$r_comp_unmgd[reef_df_indx]), c(0.025, 0.975), na.rm = TRUE)
  sector_df$ci_r_c_unmgd_lwr[sec_df_indx] <- CI[1]
  sector_df$ci_r_c_unmgd_upr[sec_df_indx] <- CI[2]
  
  # Calculate probability of recovery from compound disturbance (managed)
  ave_r_c_mgd <- reef_df$r_comp_mgd[reef_df_indx] %>%
    as.numeric() * (1 + mgmt_benefit) %>%
    min(1) 
  sector_df$ave_r_c_mgd[sec_df_indx] <- mean(ave_r_c_mgd, na.rm=TRUE)
  
  # Calculate confidence interval of recovery from single disturbance (managed)
  CI <- quantile(ave_r_c_mgd, c(0.025, 0.975), na.rm = TRUE)
  sector_df$ci_r_c_mgd_lwr[sec_df_indx] <- CI[1]
  sector_df$ci_r_c_mgd_upr[sec_df_indx] <- CI[2] # UP TO HERE
  
  # Calculate probability of recovery from single disturbance (unmanaged)
  sector_df$pr_recov_sing_unmgd[sec_df_indx] <- mean(as.numeric(reef_df$pr_recov_sing_unmgd[reef_df_indx]),
                                                     na.rm=TRUE)
  
  # Calculate probability of recovery from single disturbance (managed)
  sector_df$pr_recov_sing_mgd[sec_df_indx] <- mean(as.numeric(reef_df$pr_recov_sing_mgd[reef_df_indx]),
                                                   na.rm=TRUE)
  
  # Calculate probability of recovery from compound disturbance (unmanaged)
  sector_df$pr_recov_comp_unmgd[sec_df_indx] <- mean(as.numeric(reef_df$pr_recov_comp_unmgd[reef_df_indx]),
                                                     na.rm=TRUE)
  
  # Calculate probability of recovery from compound disturbance (managed)
  sector_df$pr_recov_comp_mgd[sec_df_indx] <- mean(as.numeric(reef_df$pr_recov_comp_mgd[reef_df_indx]),
                                                   na.rm=TRUE)
}


############################################

########## MGMT AREA OPTIMISATION ##########
# Optimise over management areas
sing_result <- optimiserSingle(sector_df, mgmt_constraint)   # Single disturbance solution
comp_result <- optimiserCompound(sector_df, mgmt_constraint) # Compound disturbance solution

# PRINT RESULTS
reef_s_to_manage_s <- which(get_solution(sing_result, y[i])$value == 1)
prcvd_exp_recov_s <- sum(sector_df$pr_recov_sing_mgd[reef_s_to_manage_s]) + 
  sum(sector_df$pr_recov_sing_unmgd[-c(reef_s_to_manage_s)])           # Perceived reefs in recov state
actual_exp_recov_s <- sum(sector_df$pr_recov_comp_mgd[reef_s_to_manage_s]) + 
  sum(sector_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_s)])           # Actual reefs in recov state
reef_s_to_manage_c <- which(get_solution(comp_result, y[i])$value == 1)
exp_recov_c <- sum(sector_df$pr_recov_comp_mgd[reef_s_to_manage_c]) + 
  sum(sector_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_c)])           # Expected reefs in recov state
print(paste("When only considering single disturbances, it is recommended to manage reef/s",
            reef_s_to_manage_s, "for a perceived", prcvd_exp_recov_s,
            "reefs in a recovered state, but an actual", actual_exp_recov_s,
            "reefs in a recovered state. In the compound disturbance scenario,",
            "it is recommended to manage reef/s", reef_s_to_manage_c,
            "for an expected", exp_recov_c, "reefs in a recovered state."))
############################################

######## MGMT AREA OPT VISUALISATION #######
# Dumbbell plot
library(ggalt)

sector_df$sector <- as.factor(sector_df$sector)

dumbbell_df <- melt(sector_df[,c("sector", 
                                 "pr_recov_comp_unmgd", 
                                 "pr_recov_sing_unmgd")], 
                    id.vars = "sector")
xy <- st_coordinates(st_centroid(sector_boundaries$geometry))
gbr_name_order <- sector_boundaries$AREA_DESCR[order(xy[,"X"], xy[,"Y"],
                                                     decreasing = TRUE)]

ggplot(data = dumbbell_df, 
       aes(x = value, y = sector)) +
  geom_line(color = "azure3") +
  geom_point(aes(color = variable), size = 3) +
  theme(legend.position = "bottom") +
  theme_light() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_color_manual(name = "", 
                     labels = c("Compound disturbance included", 
                                "Singular disturbance only"),
                     values = c("steelblue1", "steelblue4")) + 
  scale_y_discrete(limits = gbr_name_order[gbr_name_order %in% dumbbell_df$sector]) +
  scale_x_continuous(limits = c(0.5,1)) +
  labs(x = "Probability the Reef is Recovered",
       y = "GBRMPA Management Area") + 
  theme(legend.title = element_blank(),
        legend.position = c(.025, .975),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "azure3"),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(0, 5, 4, 3))
if (isTimeBased) {
  ggsave(paste0(out_path, "/Dumbbell_TimeBased", 
                recov_yrs, "yr", mgmt_benefit, "mgmt.png"), 
         plot = last_plot(), width=6, height=5)
} else {
  ggsave(paste0(out_path, "/Dumbbell_RecovBased", 
                recov_th, "th", mgmt_benefit, "mgmt.png"), 
         plot = last_plot(), width=6, height=5)
}
############################################

# To Do:
#   make raincloud vis for impact of compounding on prob recovered state

################# SAMPLING #################
num_samples <- 100
sample_reefs_df <- sampler(reef_df, sector_boundaries, num_samples, isTimeBased, 
                           recov_th, recov_yrs, cots_dist, cyc_dist, dhw_dist) %>%
  suppressMessages()

############################################

########### SAMPLE OPTIMISATION ############

sample_reefs_df <- p_calculator(sample_reefs_df, mgmt_benefit)

sing_result <- optimiserSingle(sample_reefs_df, mgmt_constraint)   # Single disturbance solution
comp_result <- optimiserCompound(sample_reefs_df, mgmt_constraint) # Compound disturbance solution

# PRINT RESULTS
reef_s_to_manage_s <- which(get_solution(sing_result, y[i])$value == 1)
prcvd_exp_recov_s <- sum(sample_reefs_df$pr_recov_sing_mgd[reef_s_to_manage_s]) + 
  sum(sample_reefs_df$pr_recov_sing_unmgd[-c(reef_s_to_manage_s)])           # Perceived reefs in recov state
actual_exp_recov_s <- sum(sample_reefs_df$pr_recov_comp_mgd[reef_s_to_manage_s]) + 
  sum(sample_reefs_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_s)])           # Actual reefs in recov state
reef_s_to_manage_c <- which(get_solution(comp_result, y[i])$value == 1)
exp_recov_c <- sum(sample_reefs_df$pr_recov_comp_mgd[reef_s_to_manage_c]) + 
  sum(sample_reefs_df$pr_recov_comp_unmgd[-c(reef_s_to_manage_c)])           # Expected reefs in recov state
print(paste("When only considering single disturbances, it is recommended to manage reef/s",
            paste(reef_s_to_manage_s, collapse = ", "), "for a perceived", 
            prcvd_exp_recov_s, "reefs in a recovered state, but an actual", 
            actual_exp_recov_s, "reefs in a recovered state.", 
            "In the compound disturbance scenario, it is recommended to manage reef/s", 
            paste(reef_s_to_manage_c, collapse = ", "), "for an expected", 
            exp_recov_c, "reefs in a recovered state."))
############################################

############### SIMULATIONS ################

# Number of reefs must be divisible by number of sections in reef_df
if (n%%length(unique(reef_df$sector)) != 0) {
  n <- floor(n/length(unique(reef_df$sector)))*
    length(unique(reef_df$sector))
}

column_names <- c('geometry', 'point_id', 'd_single', 'd_comp', 'r_single',
                  'r_comp', 'r_single_unmgd', 'r_single_mgd', 'r_comp_unmgd',
                  'r_comp_mgd', 'pr_recov_sing_unmgd', 'pr_recov_sing_mgd',
                  'pr_recov_comp_unmgd', 'pr_recov_comp_mgd')

# Initialise:
all_samples <- matrix(NA, nrow = n*n_sims, ncol = length(column_names)) %>%
  data.frame()                    # df for all sample reef info
colnames(all_samples) <- column_names
all_samples <- all_samples %>% 
  mutate(isManaged_Single = NA,   # column in df for managed/not (single)
         isManaged_Compound = NA) # column for managed/not (compound)

# For each simulation,
for (sim in 1:n_sims) {
  # Sample n reefs from reef_df
  sample_reefs_df <- sampler(reef_df, sector_boundaries, n, isTimeBased, 
                             recov_th, recov_yrs, cots_dist, cyc_dist, dhw_dist)
  
  # Calculate p for the n reefs
  sample_reefs_df <- p_calculator(sample_reefs_df, mgmt_benefit)
  
  from_row <- (sim - 1)*n + 1
  to_row <- sim*n
  all_samples[from_row:to_row, 1:length(column_names)] <- sample_reefs_df %>% 
    subset(select = c(geometry, point_id, d_single, d_comp, r_single, 
                      r_comp, r_single_unmgd, r_single_mgd, r_comp_unmgd, 
                      r_comp_mgd, pr_recov_sing_unmgd, pr_recov_sing_mgd, 
                      pr_recov_comp_unmgd, pr_recov_comp_mgd))
  
  # Calculate single disturbance solution
  sing_result <- optimiserSingle(sample_reefs_df, mgmt_constraint)   
  
  # Calculate compound disturbance solution
  comp_result <- optimiserCompound(sample_reefs_df, mgmt_constraint) 
  
  # Get the reefs to manage from solution
  sol_s <- get_solution(sing_result, y[i])$value
  sol_c <- get_solution(comp_result, y[i])$value
  
  all_samples[from_row:to_row, 'isManaged_Single'] <- sol_s
  all_samples[from_row:to_row, 'isManaged_Compound'] <- sol_c
}

# Visualise the reefs to manage
all_samples <- st_as_sf(all_samples, crs = st_crs(sector_boundaries))
single_vals <- st_transform(all_samples[all_samples$isManaged_Single == 1,], 
                            crs = st_crs(sector_boundaries))

x_single <- st_coordinates(single_vals$geometry)[,1]
y_single <- st_coordinates(single_vals$geometry)[,2]

single_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(x = "Longitude",
       y = "Latitude",
       tag = "A") +
  geom_hex(single_vals, 
           mapping = aes(x_single, y_single),
           binwidth = c(0.5,0.5)) +
  theme(plot.tag = element_text())

compound_vals <- st_transform(all_samples[all_samples$isManaged_Compound == 1,], 
                              crs = st_crs(sector_boundaries))

x_compound <- st_coordinates(compound_vals$geometry)[,1]
y_compound <- st_coordinates(compound_vals$geometry)[,2]

compound_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(x = "Longitude",
       y = "Latitude",
       tag = "B") +
  geom_hex(compound_vals, 
           mapping = aes(x_compound, y_compound),
           binwidth = c(0.5,0.5)) + 
  theme(plot.tag = element_text())
compound_plot

p <- ggplot_build(single_plot)$data[[2]]
q <- ggplot_build(compound_plot)$data[[2]]

single_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(x = "Longitude",
       y = "Latitude",
       tag = "A") +
  theme(plot.tag = element_text()) +
  geom_hex(single_vals, 
           mapping = aes(x_single, y_single),
           binwidth = c(0.5,0.5)) +
  scale_fill_continuous(limits = c(min(p$count, q$count), 
                                   max(p$count, q$count))) +
  labs(fill = "Number of \nManaged Reefs")

compound_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(x = "Longitude",
       y = "Latitude",
       tag = "B") +
  geom_hex(compound_vals, 
           mapping = aes(x_compound, y_compound),
           binwidth = c(0.5,0.5)) + 
  scale_fill_continuous(limits = c(min(p$count, q$count), 
                                   max(p$count, q$count))) +
  theme(plot.tag = element_text())

ggarrange(single_plot, compound_plot, 
          ncol=2, nrow=1, 
          common.legend = TRUE, 
          legend="right")

if (isTimeBased) {
  ggsave(paste0(out_path, "/", 
                n_sims, "Sim_Hexmap_TimeBased", 
                recov_yrs, "yr", mgmt_benefit, "mgmt.png"), 
         plot = last_plot(), 
         width=8, height=5)
} else {
  ggsave(paste0(out_path, "/", 
                n_sims, "Sim_Hexmap_RecovBased", 
                recov_th, "th", mgmt_benefit, "mgmt.png"), 
         plot = last_plot(), 
         width=8, height=5)
}


############################################

############# SAVE ENVIRONMENT #############
if (isTimeBased) {
  save.image(paste0(data_path,
                    "/Parameter Run Environments/", 
                    "TimeBased", recov_yrs, "yrs", 
                    mgmt_benefit, "mgmt.RData"))
} else {
  save.image(paste0(data_path,
                    "/Parameter Run Environments/", 
                    "RecovBased", recov_th, "th", 
                    mgmt_benefit, "mgmt.RData"))
}
############################################

############ LIBRARY GRAVEYARD #############
# library(ncdf4)            # opens ".nc" files
# library(colorspace)       # makes colouring plots easier
# library(lubridate)        # manipulates date/time easily
# library(tidyverse)        # pipes
# library(stringr)          # handles strings
# library(lpSolve)          # linear programming model solver
# library(leaflet)
# library(maps)
# library(raster)
# library(ggspatial)
# library(dataaimsr)
# library(gisaimsr)
# library(ggrepel)
# library(leaflet.providers)
############################################

############# CODE GRAVEYARD ###############

############################################

############### VIS GRAVEYARD ##############

# SAMPLE OPT VIS
# Plot the sectors with different colours
# sector_boundaries$AREA_DESCR <- factor(sector_boundaries$AREA_DESCR, 
#                                        levels = gbr_name_order[4:1])
# 
# sample_reefs_df <- st_as_sf(sample_reefs_df, crs = st_crs(sector_boundaries))
# 
# sample_reefs_df$x <- st_coordinates(sample_reefs_df$geometry)[,1]
# sample_reefs_df$y <- st_coordinates(sample_reefs_df$geometry)[,2]
# 
# sector_boundaries_fill <- sector_boundaries %>% 
#   mutate(single_count = nrow(sample_reefs_df[sample_reefs_df$AREA_DESCR[reef_s_to_manage_s] == 
#                                                AREA_DESCR,]), 
#          compound_count = nrow(sample_reefs_df[sample_reefs_df$AREA_DESCR[reef_s_to_manage_c] == 
#                                                  AREA_DESCR,]))
# for (sector in sector_boundaries_fill$AREA_DESCR) {
#   sector_boundaries_fill$single_count[sector_boundaries_fill$AREA_DESCR == sector] = 
#     sum(sample_reefs_df$AREA_DESCR[reef_s_to_manage_s] == sector)
#   sector_boundaries_fill$compound_count[sector_boundaries_fill$AREA_DESCR == sector] = 
#     sum(sample_reefs_df$AREA_DESCR[reef_s_to_manage_c] == sector)
# }
# 
# single_p <- ggplot() +
#   geom_sf(data = sector_boundaries_fill,
#           mapping = aes(fill = single_count), lwd = 0.01) +
#   theme_classic() +
#   labs(x = "Longitude",
#        y = "Latitude",
#        tag = "A") +
#   scale_fill_continuous(limits = c(min(sector_boundaries_fill$single_count, 
#                                        sector_boundaries_fill$compound_count), 
#                                    max(sector_boundaries_fill$single_count, 
#                                        sector_boundaries_fill$compound_count))) +
#   labs(fill = "Number of \nManaged Reefs") +
#   theme(plot.tag = element_text()) + 
#   annotate("text", x = 145, y = -24, 
#            label = TeX(paste0("$E\\[R_{1}\\] = ", floor(actual_exp_recov_s), "$")),
#            parse = TRUE) %>%
#   suppressMessages()
# 
# compound_p <- ggplot() +
#   geom_sf(data = sector_boundaries_fill,
#           mapping = aes(fill = compound_count), lwd = 0.01) +
#   theme_classic() +
#   labs(x = "Longitude",
#        y = "Latitude",
#        tag = "B") +
#   scale_fill_continuous(limits = c(min(sector_boundaries_fill$single_count, 
#                                        sector_boundaries_fill$compound_count), 
#                                    max(sector_boundaries_fill$single_count, 
#                                        sector_boundaries_fill$compound_count))) +
#   theme(plot.tag = element_text()) +
#   annotate("text", x = 145, y = -24, 
#            label = TeX(paste0("$E\\[R_{2}\\] = ", floor(exp_recov_c), "$")),
#            parse = TRUE) %>%
#   suppressMessages()
# 
# ggarrange(single_p, compound_p, 
#           ncol=2, nrow=1, 
#           common.legend = TRUE, 
#           legend="right") %>%
#   suppressMessages()
# 
# if (isTimeBased) {
#   ggsave(paste0("../MPhil Thesis/", 
#                 num_samples, "Sample_Choropleth_TimeBased", 
#                 recov_yrs, "yr", mgmt_benefit, "mgmt.png"), 
#          plot = last_plot(), 
#          width=8, height=5)
# } else {
#   ggsave(paste0("../MPhil Thesis/", 
#                 num_samples, "Sample_Choropleth_RecovBased", 
#                 recov_th, "th", mgmt_benefit, "mgmt.png"), 
#          plot = last_plot(), 
#          width=8, height=5)
# }

# # Locations of (real) reefs in system
# ggplot() +
#   geom_sf(data = sector_boundaries,
#           mapping = aes(fill = AREA_DESCR), lwd = 0.01) +
#   theme_classic() +
#   labs(x = "Longitude",
#        y = "Latitude") +
#   scale_fill_brewer(name = "Region", palette = "Spectral") +
#   geom_sf(data = all_reefs_sf$geometry,
#           pch = 16,
#           color = "black",
#           alpha = 0.5)
############################################