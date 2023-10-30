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

############################################

############## SET VARIABLES ###############
# Variables
base_recov_th <- 0.75
base_mgmt_ben <- 0.5
list_recov_th <- seq(0, 1, by = 0.1)
list_mgmt_ben <- seq(0, 1, by = 0.1)
list_recov_yrs <- seq(1, 15, by = 1)
infer_baseline <- 1 # Should the baseline be inferred if non-existent at start of obs period?
epsilon <- 0.05
baseline_str <- "mean"

# Management constraint (base: 20% of number of reefs in system)
mgmt_constraint <- 0.20

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
out_path <- paste0(mphil_path, "/Figures/OptChapter/Sensitivity")

# Time based or recovery based compounding?
#  Note: Set to TRUE for time-based compounding or FALSE for recovery-based
is_time_based <- FALSE

### Simulations ###
# Run simulations? (Much faster if you don't)
run_simulations <- TRUE

# Set number of simulations, n_sims
n_sims <- 500

# Set number of sample reefs, num_samples
num_samples <- 100

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

# SIMULATIONS
column_names <- c(
    "sector", "point_loc", "point_id", "num_dist", "num_s_dist", "num_c_dist",
    "prob_s_dist", "prob_c_dist", "prob_s_impact", "prob_c_impact",
    "prob_s_recov", "prob_c_recov", "r_single_unmgd", "r_single_mgd",
    "r_comp_unmgd", "r_comp_mgd", "pr_recov_sing_unmgd", "pr_recov_sing_mgd",
    "pr_recov_comp_unmgd", "pr_recov_comp_mgd"
)

# Make a dataframe to store the results
opt_sensitivity_1 <- expand.grid(
    recov_th = list_recov_th,
    mgmt_benefit = list_mgmt_ben
) %>% mutate(
    expected_recov_s = NA,
    expected_recov_c = NA,
    expected_recov_nm = NA
)

# Make a list to store all simulations
sens_list_1 = replicate(n = nrow(opt_sensitivity_1),
                        expr = {matrix(nrow = num_samples * n_sims,
                                       ncol = length(column_names),
                                       dimnames = list(NULL, column_names)) %>%
                                data.frame() %>%
                                mutate(recov_th = NA,
                                       mgmt_benefit = NA)
                        },
                        simplify = F)


############################################

########### SENSITIVITY 1 LOOP #############
i = 1
total = length(list_recov_th) * length(list_mgmt_ben)
for (recov_th in list_recov_th) {
    for (mgmt_benefit in list_mgmt_ben) {
        print(paste("Starting i =", i, "out of 121",
                    "recov_th =", recov_th, 
                    ", mgmt_benefit =", mgmt_benefit))
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
        reef_df <- sing_comp_list[[3]]
        
        # Calculate probability of reef being in a recovered state
        reef_df <- p_calculator(reef_df, mgmt_benefit)
        
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
            
            if(any(!is.na(sample_reefs_df$pr_recov_sing_unmgd) &
                   !is.na(sample_reefs_df$pr_recov_sing_mgd))) {
                # Calculate single disturbance solution
                sing_result <- optimiser_single(sample_reefs_df, mgmt_constraint)
                
                # Get the reefs to manage from solution
                sol_s <- get_solution(sing_result, y[i])$value
                all_samples[from_row:to_row, "is_managed_single"] <- sol_s
            }
            
            
            if(any(!is.na(sample_reefs_df$pr_recov_comp_unmgd) &
                   !is.na(sample_reefs_df$pr_recov_comp_mgd))) {
                # Calculate compound disturbance solution
                comp_result <- optimiser_compound(sample_reefs_df, mgmt_constraint)
                
                # Get the reefs to manage from solution
                sol_c <- get_solution(comp_result, y[i])$value
                all_samples[from_row:to_row, "is_managed_cumul"] <- sol_c
            }
        }
        
        # Save all_samples in list
        sens_list_1[[i]] <- all_samples %>%
            mutate(recov_th = recov_th,
                   mgmt_benefit = mgmt_benefit)
        i <- i + 1
        
        # Calculate the expected number of reefs recovered
        expected_recov_reefs_s <- all_samples %>%
            group_by(sim_num) %>%
            summarise(num_recov_s = 
                          sum(pr_recov_comp_mgd[all_samples$is_managed_single == 1],
                              na.rm = TRUE) +
                          sum(pr_recov_comp_unmgd[all_samples$is_managed_single == 0],
                              na.rm = TRUE)
            )
        
        expected_recov_reefs_c <- all_samples %>%
            group_by(sim_num) %>%
            summarise(num_recov_c = 
                          sum(pr_recov_comp_mgd[all_samples$is_managed_cumul == 1],
                              na.rm = TRUE) +
                          sum(pr_recov_comp_unmgd[all_samples$is_managed_cumul == 0],
                              na.rm = TRUE)
            )
        
        expected_recov_reefs_nm <- all_samples %>%
            group_by(sim_num) %>%
            summarise(num_recov_nm = sum(pr_recov_comp_unmgd,
                                         na.rm = TRUE)
            )
        
        # Save the results
        sensitivity_row <- which(opt_sensitivity_1$recov_th == recov_th &
                                     opt_sensitivity_1$mgmt_benefit == mgmt_benefit)
        opt_sensitivity_1[sensitivity_row, "expected_recov_s"] <- mean(expected_recov_reefs_s$num_recov_s)
        opt_sensitivity_1[sensitivity_row, "expected_recov_c"] <- mean(expected_recov_reefs_c$num_recov_c)
        opt_sensitivity_1[sensitivity_row, "expected_recov_nm"] <- mean(expected_recov_reefs_nm$num_recov_nm)
    }
}
############################################

########### FIG 1: OPT TO PARAMS ###########
opt_sensitivity_1 <- opt_sensitivity_1 %>%
    mutate(difference_s = expected_recov_s - expected_recov_nm,
           difference_c = expected_recov_c - expected_recov_nm,
           prop_benefit = ifelse(is.nan(difference_c / difference_s),
                                 0,
                                 difference_c / difference_s))

ggplot(opt_sensitivity_1, aes(x = as.factor(recov_th*100), 
                              y = as.factor(mgmt_benefit*100), 
                              fill = (prop_benefit))) + 
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white",
                         high = "blue", space = "Lab" ) +
    labs(fill = "Propotional benefit of incorporating\nsingle and cumulative disturbance\nfor expected number of recovered reefs\n(multiple of general benefit)",
         x = "Recovery Threshold (%)",
         y = "Management Benefit (%)") +
    theme_classic()

ggsave(paste0(out_path, "/Sensitivity1.png"), 
       plot = last_plot(), 
       width=8, height=5)

############################################

############# SAVE ENVIRONMENT #############
save.image(paste0(
    data_path,
    "/SensitivityData/Opt/Sens1.RData"
))
############################################

########### SENSITIVITY 2 LOOP #############
bc_opt_indx <- which(round(opt_sensitivity_1$recov_th, 2) == base_recov_th &
                         opt_sensitivity_1$mgmt_benefit == base_mgmt_ben)
base_case_optimal <- opt_sensitivity_1[bc_opt_indx, "expected_recov_c"]
is_time_based <- TRUE
mgmt_benefit <- base_mgmt_ben

# Load disturbances and coral obs (.rds made in code/obs_dist_combine.R)
all_reefs_sf <- readRDS(file = paste0(
    data_path,
    "/all_reefs_sf_gaps_filled.rds"
))

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

# Make a dataframe to store the results
opt_sensitivity_2 <- data.frame(
    recov_yrs = list_recov_yrs,
    expected_recov_s = NA,
    expected_recov_c = NA
)

# Make a list to store all simulations
sens_list_2 = replicate(n = nrow(opt_sensitivity_2),
                        expr = {matrix(nrow = num_samples * n_sims,
                                       ncol = length(column_names),
                                       dimnames = list(NULL, column_names)) %>%
                                data.frame() %>%
                                mutate(recov_yrs = NA)
                        },
                        simplify = F)
i = 1
for (recov_yrs in list_recov_yrs) {
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
    reef_df <- sing_comp_list[[3]]
    
    # Calculate probability of reef being in a recovered state
    reef_df <- p_calculator(reef_df, mgmt_benefit)
    
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
    
    # Save all_samples in list
    sens_list_2[[i]] <- all_samples %>%
        mutate(recov_yrs = recov_yrs)
    i <- i + 1
    
    # Calculate the expected number of reefs recovered
    expected_recov_reefs_s <- all_samples %>%
        group_by(sim_num) %>%
        summarise(num_recov_s = 
                      sum(pr_recov_comp_mgd[all_samples$is_managed_single == 1],
                          na.rm = TRUE) +
                      sum(pr_recov_comp_unmgd[all_samples$is_managed_single == 0],
                          na.rm = TRUE)
        )
    
    expected_recov_reefs_c <- all_samples %>%
        group_by(sim_num) %>%
        summarise(num_recov_c = 
                      sum(pr_recov_comp_mgd[all_samples$is_managed_cumul == 1],
                          na.rm = TRUE) +
                      sum(pr_recov_comp_unmgd[all_samples$is_managed_cumul == 0],
                          na.rm = TRUE)
        )
    
    # Save the results
    sensitivity_row <- which(opt_sensitivity_2$recov_yrs == recov_yrs)
    opt_sensitivity_2[sensitivity_row, "expected_recov_s"] <- mean(expected_recov_reefs_s$num_recov_s)
    opt_sensitivity_2[sensitivity_row, "expected_recov_c"] <- mean(expected_recov_reefs_c$num_recov_c)
}

############################################


########## FIG 2: OPT TO TIME EST ##########
opt_sensitivity_2 <- opt_sensitivity_2 %>%
    mutate(difference = abs(base_case_optimal - expected_recov_c))

ggplot(opt_sensitivity_2, aes(x = as.factor(recov_yrs), 
                              y = difference)) + 
    geom_point() +
    labs(x = "Estimated Recovery Time Following Disturbance (Years)",
         y = "Difference in average expected number of recovered reefs") +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 45, 5), limits = c(0, 45))

ggsave(paste0(out_path, "/Sensitivity2.png"), 
       plot = last_plot(), 
       width=5, height=5)

# opt_sensitivity_2 <- sens_list_2 %>%
#     bind_rows()
# 
# # Make a new column for the probability of recovery based on whether the reef is managed or not
# opt_sensitivity_2 <- opt_sensitivity_2 %>%
#     mutate(pr_recov_sing = ifelse(is_managed_single == 1, pr_recov_comp_mgd, pr_recov_comp_unmgd),
#            pr_recov_cumul = ifelse(is_managed_cumul == 1,  pr_recov_comp_mgd, pr_recov_comp_unmgd))
# opt_sensitivity_2_fig <- opt_sensitivity_2 %>%
#     group_by(sim_num, recov_yrs) %>%
#     summarise(num_recov_s = 
#                   sum(pr_recov_sing,
#                       na.rm = TRUE),
#               num_recov_c =
#                   sum(pr_recov_cumul,
#                       na.rm = TRUE)
#     )
# 
# ggplot(opt_sensitivity_2_fig, aes(x = as.factor(recov_yrs), 
#                                   y = num_recov_c - num_recov_s)) + 
#     geom_boxplot() +
#     labs(x = "Estimated Recovery Time Following Disturbance (Years)",
#          y = "Difference in average expected number of recovered reefs") +
#     theme_classic()
# 
# ggsave(paste0(out_path, "/Sensitivity2.png"), 
#        plot = last_plot(), 
#        width=5, height=5)

############################################

############# SAVE ENVIRONMENT #############
save.image(paste0(
    data_path,
    "/SensitivityData/Opt2.RData"
))
############################################
