########### SET UP THE WORKSPACE ###########
# Clear environment
rm(list = ls())

# Clear plots
if (!is.null(dev.list())) dev.off()

# Clear commands
cat("\014")

############## LOAD LIBRARIES ###############
# Check if the required packages are installed and load them
required_packages <- c("sf", "lwgeom", "dplyr", "varhandle", "ompr", 
                       "ompr.roi", "reshape2", "nngeo", "lubridate", 
                       "tidyr", "ggpubr", "ROI.plugin.glpk", "ROI", "ggplot2", 
                       "latex2exp", "patchwork", "ggbump", "GGally", "pracma",
                       "parallel", "progress")
new_packages <- required_packages[!(required_packages %in% 
                                      installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
# Load Libraries
suppressMessages(
  suppressWarnings(
    lapply(required_packages, require, character.only = TRUE)
  )
)

# Load custom functions
source_files <- c("Stewart_MPhil_single_or_compound.R", 
                  "Stewart_MPhil_ortiz_r_func.R", 
                  "Stewart_MPhil_p_calc.R", 
                  "Stewart_MPhil_optimiser_single.R", 
                  "Stewart_MPhil_optimiser_compound.R", 
                  "sampler_v2.R", 
                  "Stewart_MPhil_dist_finder.R", 
                  "Stewart_MPhil_analyser.R")

sapply(source_files, source)

# Set paths
qut_path <- "../OneDrive - Queensland University of Technology"
mphil_path <- paste0(qut_path, "/Documents/MPhil")
data_path <- paste0(mphil_path, "/Data")
out_path <- paste0(mphil_path, "/Figures/OptChapter/Sensitivity")

# Load data
all_reefs_sf <- readRDS(file = paste0(data_path, "/all_reefs_sf_gaps_filled.rds"))
sample_reefs <- readRDS(file = paste0(data_path, "/sample_reefs.rds"))
shapefile_path <- paste0(data_path, "/GBRMPA/Management_Areas_of_the_Great_Barrier",
                         "_Reef_Marine_Park/Management_Areas_of_the_Great_",
                         "Barrier_Reef_Marine_Park.shp")
sector_boundaries <- st_read(shapefile_path, quiet = TRUE)

############## SET VARIABLES ###############
base_recov_th <- 0.25
base_mgmt_ben <- 0.5
list_recov_th <- seq(0, 1, by = 0.05)
list_mgmt_ben <- seq(0, 1, by = 0.05)
list_recov_yrs <- seq(1, 15, by = 5)
infer_baseline <- 1
epsilon <- 0.05
baseline_str <- "max"
mgmt_constraint <- 0.20

manament_area_names <- c(
  "Far Northern Management Area",
  "Cairns/Cooktown Management Area",
  "Townsville/Whitsunday Management Area",
  "Mackay/Capricorn Management Area"
)

is_time_based <- FALSE
run_simulations <- TRUE
n_sims <- 1000
num_samples <- 100
cots_dist <- 1
cyc_dist <- 10
dhw_dist <- 4

reef_names <- unique(all_reefs_sf$REEF_NAME)
overall_unknown_count <- 0

all_reefs_sf <- all_reefs_sf %>%
  mutate(
    is_disturbed = (COTS_value >= cots_dist |
                      Hs4MW_value >= cyc_dist |
                      DHW_value >= dhw_dist),
    is_impacted = 0,
    single_or_compound = NA,
    dist_type = NA,
    recov_year = NA,
    recov_time = NA,
    r_given_impact = NA
  )

# Initialize results dataframes
opt_sensitivity_1 <- expand.grid(
  recov_th = list_recov_th,
  mgmt_benefit = list_mgmt_ben
) %>%
  mutate(
    expected_recov_s = NA,
    expected_recov_c = NA,
    expected_recov_nm = NA
  )

opt_sensitivity_2 <- data.frame(
  recov_yrs = list_recov_yrs,
  expected_recov_s = NA,
  expected_recov_c = NA
)

column_names <- c(
  "sector", "point_loc", "point_id", "num_dist", "num_s_dist", "num_c_dist",
  "prob_s_dist", "prob_c_dist", "prob_s_impact", "prob_c_impact",
  "prob_s_recov", "prob_c_recov", "r_single_unmgd", "r_single_mgd",
  "r_comp_unmgd", "r_comp_mgd", "pr_recov_sing_unmgd", "pr_recov_sing_mgd",
  "pr_recov_comp_unmgd", "pr_recov_comp_mgd"
)

# Define utility function
run_analysis <- function(recov_th, mgmt_benefit, all_reefs_sf, num_samples, 
                         n_sims, column_names, mgmt_constraint, 
                         sector_boundaries, sample_reefs, is_time_based, 
                         recov_yrs, cots_dist, cyc_dist, dhw_dist, 
                         infer_baseline, epsilon, baseline_str) {
  sing_comp_list <- analyse_reefs(
    all_reefs_sf, reef_names,
    is_time_based, recov_yrs, recov_th, cots_dist,
    cyc_dist, dhw_dist, infer_baseline, epsilon,
    baseline_str
  )
  
  event_counts <- sing_comp_list[[1]]
  all_reefs_sf <- sing_comp_list[[2]]
  reef_df <- sing_comp_list[[3]]
  
  reef_df <- p_calculator(reef_df, mgmt_benefit)
  
  all_samples <- matrix(NA, 
                        nrow = num_samples * n_sims, 
                        ncol = length(column_names)) %>%
    data.frame() 
  colnames(all_samples) <- column_names
  all_samples <- all_samples %>%
    mutate(
      is_managed_single = NA,
      is_managed_cumul = NA,
      scenario_managed = NA,
      sim_num = rep(seq(n_sims), each = num_samples)
    )
  
  for (sim in 1:n_sims) {
    sample_reefs_df <- samplerv2(
      reef_df, sector_boundaries, sample_reefs, num_samples,
      is_time_based, recov_th, recov_yrs, cots_dist, cyc_dist,
      dhw_dist
    )
    
    sample_reefs_df <- p_calculator(sample_reefs_df, mgmt_benefit)
    
    from_row <- (sim - 1) * num_samples + 1
    to_row <- sim * num_samples
    all_samples[from_row:to_row, 1:length(column_names)] <- sample_reefs_df %>%
      subset(select = column_names)
    
    if(any(!is.na(sample_reefs_df$pr_recov_sing_unmgd) &
           !is.na(sample_reefs_df$pr_recov_sing_mgd))) {
      sing_result <- optimiser_single(sample_reefs_df, mgmt_constraint)
      sol_s <- get_solution(sing_result, y[i])$value
      all_samples[from_row:to_row, "is_managed_single"] <- sol_s
    }
    
    if(any(!is.na(sample_reefs_df$pr_recov_comp_unmgd) &
           !is.na(sample_reefs_df$pr_recov_comp_mgd))) {
      comp_result <- optimiser_compound(sample_reefs_df, mgmt_constraint)
      sol_c <- get_solution(comp_result, y[i])$value
      all_samples[from_row:to_row, "is_managed_cumul"] <- sol_c
    }
  }
  
  list(
    all_samples = all_samples,
    expected_recov_s = all_samples %>%
      group_by(sim_num) %>%
      summarise(num_recov_s = 
                  sum(pr_recov_comp_mgd[all_samples$is_managed_single == 1], 
                      na.rm = TRUE) +
                  sum(pr_recov_comp_unmgd[all_samples$is_managed_single == 0], 
                      na.rm = TRUE)
      ) %>% 
      summarise(mean(num_recov_s)) %>% 
      pull(),
    expected_recov_c = all_samples %>%
      group_by(sim_num) %>%
      summarise(num_recov_c = 
                  sum(pr_recov_comp_mgd[all_samples$is_managed_cumul == 1], 
                      na.rm = TRUE) +
                  sum(pr_recov_comp_unmgd[all_samples$is_managed_cumul == 0], 
                      na.rm = TRUE)
      ) %>% 
      summarise(mean(num_recov_c)) %>% 
      pull(),
    expected_recov_nm = all_samples %>%
      group_by(sim_num) %>%
      summarise(num_recov_nm = sum(pr_recov_comp_unmgd, na.rm = TRUE)
      ) %>% 
      summarise(mean(num_recov_nm)) %>% 
      pull()
  )
}

# Run sensitivity analysis in parallel
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(dplyr)
  library(sf)
  library(lwgeom)
  library(varhandle)
  library(ompr)
  library(ompr.roi)
  library(reshape2)
  library(nngeo)
  library(lubridate)
  library(tidyr)
  library(ggpubr)
  library(ROI.plugin.glpk)
  library(ROI)
  library(ggplot2)
  library(latex2exp)
  library(patchwork)
  library(ggbump)
  library(GGally)
  library(pracma)
})

# Initialize progress bar
pb <- progress_bar$new(
  format = "  Running [:bar] :percent in :elapsed",
  total = nrow(opt_sensitivity_1), clear = FALSE, width = 60
)

sensitivity_results <- parLapply(cl, 1:nrow(opt_sensitivity_1), function(i) {
  recov_th <- opt_sensitivity_1$recov_th[i]
  mgmt_benefit <- opt_sensitivity_1$mgmt_benefit[i]
  run_analysis(recov_th, mgmt_benefit, all_reefs_sf, num_samples, n_sims, 
               column_names, mgmt_constraint, sector_boundaries, sample_reefs, 
               is_time_based, recov_yrs, cots_dist, cyc_dist, dhw_dist, 
               infer_baseline, epsilon, baseline_str)
  pb$tick()
  return(result)
})

stopCluster(cl)

# Collect results
for (i in 1:nrow(opt_sensitivity_1)) {
  opt_sensitivity_1$expected_recov_s[i] <- sensitivity_results[[i]]$expected_recov_s
  opt_sensitivity_1$expected_recov_c[i] <- sensitivity_results[[i]]$expected_recov_c
  opt_sensitivity_1$expected_recov_nm[i] <- sensitivity_results[[i]]$expected_recov_nm
}

# Plot results
opt_sensitivity_1 <- opt_sensitivity_1 %>%
  mutate(difference_s = expected_recov_s - expected_recov_nm,
         difference_c = expected_recov_c - expected_recov_nm,
         prop_benefit = ifelse(is.nan(difference_c / difference_s), 0, difference_c / difference_s))

sens_1_plot <- ggplot(opt_sensitivity_1, 
                      aes(x = as.factor(recov_th*100), 
                          y = as.factor(mgmt_benefit*100), 
                          fill = (prop_benefit))) + 
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white",
                       high = "blue", space = "Lab") +
  labs(fill = "Proportional benefit of incorporating\nsingle and cumulative disturbance\nfor expected number of recovered reefs\n(multiple of general benefit)",
       x = "Recovery Threshold (%)",
       y = "Management Benefit (%)") +
  theme_classic()

ggsave(paste0(out_path, "/Sensitivity1.png"), 
       plot = sens_1_plot, 
       width=8, height=5)
############################################

############# SAVE ENVIRONMENT #############
save.image("SensitivityData/Opt1.RData")
############################################

# Sensitivity 2 Loop
bc_opt_indx <- which(round(opt_sensitivity_1$recov_th, 2) == base_recov_th &
                       opt_sensitivity_1$mgmt_benefit == base_mgmt_ben)
base_case_optimal <- opt_sensitivity_1[bc_opt_indx, "expected_recov_c"]
is_time_based <- TRUE
mgmt_benefit <- base_mgmt_ben

all_reefs_sf <- all_reefs_sf %>%
  mutate(
    is_disturbed = (COTS_value >= cots_dist |
                      Hs4MW_value >= cyc_dist |
                      DHW_value >= dhw_dist),
    is_impacted = 0,
    single_or_compound = NA,
    dist_type = NA,
    recov_year = NA,
    recov_time = NA,
    r_given_impact = NA
  )

# Run sensitivity analysis in parallel
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, {
  library(dplyr)
  library(sf)
  library(lwgeom)
  library(varhandle)
  library(ompr)
  library(ompr.roi)
  library(reshape2)
  library(nngeo)
  library(lubridate)
  library(tidyr)
  library(ggpubr)
  library(ROI.plugin.glpk)
  library(ROI)
  library(ggplot2)
  library(latex2exp)
  library(patchwork)
  library(ggbump)
  library(GGally)
  library(pracma)
})

# Run second sensitivity analysis
pb2 <- progress_bar$new(
  format = "  Running [:bar] :percent in :elapsed",
  total = nrow(opt_sensitivity_1), clear = FALSE, width = 60
)
sensitivity_results_2 <- parLapply(cl, 1:nrow(opt_sensitivity_2), function(i) {
  recov_yrs <- opt_sensitivity_2$recov_yrs[i]
  run_analysis(recov_th, mgmt_benefit, all_reefs_sf, num_samples, n_sims, 
               column_names, mgmt_constraint, sector_boundaries, sample_reefs, 
               is_time_based, recov_yrs, cots_dist, cyc_dist, dhw_dist, 
               infer_baseline, epsilon, baseline_str)
  pb2$tick()
  return(result)
})

stopCluster(cl)

# Collect results for second sensitivity analysis
for (i in 1:nrow(opt_sensitivity_2)) {
  opt_sensitivity_2$expected_recov_s[i] <- sensitivity_results_2[[i]]$expected_recov_s
  opt_sensitivity_2$expected_recov_c[i] <- sensitivity_results_2[[i]]$expected_recov_c
}

# Plot second sensitivity analysis results
opt_sensitivity_2 <- opt_sensitivity_2 %>%
  mutate(difference = abs(base_case_optimal - expected_recov_c))

sens_2_plot <- ggplot(opt_sensitivity_2, aes(x = as.factor(recov_yrs), y = difference)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(x = "Estimated Recovery Time Following Disturbance (Years)",
       y = "Difference in average expected number of recovered reefs") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 45, 5), limits = c(0, 45))

ggsave(paste0(out_path, "/Sensitivity2.png"), 
       plot = sens_2_plot, 
       width=5, height=5)

############# SAVE ENVIRONMENT #############
save.image("SensitivityData/Opt2.RData")
############################################
