analyse_reefs <- function(
        all_reefs_sf, reef_names, is_time_based,
        recov_yrs, recov_th, cots_dist, cyc_dist,
        dhw_dist, infer_baseline, epsilon, baseline_str) {
    # Initialise dataframe of disturbance type and counts
    event_counts <- data.frame(
        H = 0,
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
        Other = 0
    )
    
    # Initialise dataframe for reef information
    reef_df <- data.frame(
        reef_name = character(),
        latitude = double(),
        longitude = double(),
        sector = character(),
        num_total = integer(),
        num_single = integer(),
        prob_s_dist = double(),
        prob_s_impact = double(),
        prob_s_recov = double(),
        num_comp = integer(),
        prob_c_dist = double(),
        prob_c_impact = double(),
        prob_c_recov = double(),
        reef_unknown = integer(),
        baseline_inferred = character(),
        baseline_vals = character()
    )
    
    # Initialise overall unknown count
    overall_unknown_count <- 0
    
    # Run code for each reef
    for (reef_name in reef_names) {
        # Extract reef obs for named reef
        obs_by_reef <- all_reefs_sf[all_reefs_sf$REEF_NAME == reef_name, ]
        # Separate single and compound disturbances and get recovery times
        if (is_time_based) {
            result_vec <- single_or_compound(
                obs_by_reef = obs_by_reef,
                is_time_based = is_time_based,
                recov_yrs = recov_yrs,
                cots_dist = cots_dist,
                cyc_dist = cyc_dist,
                dhw_dist = dhw_dist,
                infer_baseline = infer_baseline,
                epsilon = epsilon,
                baseline_str = baseline_str
            )
            
            obs_by_reef <- result_vec[[1]]
        } else {
            result_vec <- single_or_compound(
                obs_by_reef = obs_by_reef,
                is_time_based = is_time_based,
                recov_th = recov_th,
                cots_dist = cots_dist,
                cyc_dist = cyc_dist,
                dhw_dist = dhw_dist,
                infer_baseline = infer_baseline,
                epsilon = epsilon,
                baseline_str = baseline_str
            )
            
            obs_by_reef <- result_vec[[1]]
        }
        
        if (!is.null(result_vec[[2]])) {
            baseline_inferred <- result_vec[[2]]
        } else {
            baseline_inferred <- FALSE
        }
        if (!is.null(result_vec[[3]])) {
            baseline_vals <- result_vec[[3]]
        } else {
            baseline_vals <- NULL
        }
        
        # Calculate number of years observed inclusive
        yrs_obsvd <- obs_by_reef$YEAR[nrow(obs_by_reef)] - obs_by_reef$YEAR[1] + 1
        # Calculate disturbance probabilities given occurrences
        # d = (number of events) / (years observed)
        single_dist <- obs_by_reef %>%
            filter(single_or_compound == "Single")
        num_single <- nrow(single_dist)
        prob_s_dist <- num_single / yrs_obsvd
        # Prob Impacted given disturbance: P(A|B) = P(A & B)/P(B)
        prob_s_impact <- ifelse(num_single == 0,
                                NA,
                                sum(single_dist$is_impacted & single_dist$is_disturbed,
                                    na.rm = TRUE
                                ) /
                                    sum(single_dist$is_disturbed,
                                        na.rm = TRUE
                                    )
        )
        if (!is.na(prob_s_impact) && prob_s_impact < 0) {
            prob_s_impact <- 0
        } else if (!is.na(prob_s_impact) && prob_s_impact > 1) {
            prob_s_impact <- 1
        }
        prob_s_recov <- ifelse(any(!is.na(single_dist$r_given_impact)) && 
                                   !is.na(sd(single_dist$r_given_impact, 
                                             na.rm = TRUE)),
                               rnorm(1, 
                                     mean(single_dist$r_given_impact, 
                                          na.rm = TRUE), 
                                     sd(single_dist$r_given_impact, 
                                        na.rm = TRUE)
                               ),
                               NA
        )
        if (!is.na(prob_s_recov) && prob_s_recov < 0) {
            prob_s_recov <- 0
        } else if (!is.na(prob_s_recov) && prob_s_recov > 1) {
            prob_s_recov <- 1
        }
        comp_dist <- obs_by_reef %>%
            filter(single_or_compound == "Compound")
        num_comp <- nrow(comp_dist)
        prob_c_dist <- num_comp / yrs_obsvd
        prob_c_impact <- ifelse(num_comp == 0,
                                NA,
                                sum(comp_dist$is_impacted & comp_dist$is_disturbed,
                                    na.rm = TRUE
                                ) /
                                    sum(comp_dist$is_disturbed,
                                        na.rm = TRUE
                                    )
        )
        if (!is.na(prob_c_impact) && prob_c_impact < 0) {
            prob_c_impact <- 0
        } else if (!is.na(prob_c_impact) && prob_c_impact > 1) {
            prob_c_impact <- 1
        }
        prob_c_recov <- ifelse(any(!is.na(comp_dist$r_given_impact)) && 
                                   !is.na(sd(comp_dist$r_given_impact[!is.na(comp_dist$r_given_impact)], 
                                             na.rm = TRUE)),
                               rnorm(1, 
                                     mean(comp_dist$r_given_impact[!is.na(comp_dist$r_given_impact)], na.rm = TRUE), 
                                     sd(comp_dist$r_given_impact[!is.na(comp_dist$r_given_impact)], na.rm = TRUE)
                               ),
                               NA
        )
        if (!is.na(prob_c_recov) && prob_c_recov < 0) {
            prob_c_recov <- 0
        } else if (!is.na(prob_c_recov) && prob_c_recov > 1) {
            prob_c_recov <- 1
        }
        # Get count of events for the reef and add to total event_counts df
        event_counts_rf <- c(
            H = sum(obs_by_reef$dist_type == "Heat Stress",
                    na.rm = TRUE
            ),
            Cy = sum(obs_by_reef$dist_type == "Wind Stress",
                     na.rm = TRUE
            ),
            Co = sum(obs_by_reef$dist_type == "CoTS Outbreak",
                     na.rm = TRUE
            ),
            HH = sum(obs_by_reef$dist_type == "Heat Stress, Heat Stress",
                     na.rm = TRUE
            ),
            HCy = sum(obs_by_reef$dist_type == "Heat Stress, Wind Stress",
                      na.rm = TRUE
            ),
            CyH = sum(obs_by_reef$dist_type == "Wind Stress, Heat Stress",
                      na.rm = TRUE
            ),
            CyCy = sum(obs_by_reef$dist_type == "Wind Stress, Wind Stress",
                       na.rm = TRUE
            ),
            CyCo = sum(obs_by_reef$dist_type == "Wind Stress, CoTS Outbreak",
                       na.rm = TRUE
            ),
            CoCy = sum(obs_by_reef$dist_type == "CoTS Outbreak, Wind Stress",
                       na.rm = TRUE
            ),
            CoCo = sum(obs_by_reef$dist_type == "CoTS Outbreak, CoTS Outbreak",
                       na.rm = TRUE
            ),
            CoH = sum(obs_by_reef$dist_type == "CoTS Outbreak, Heat Stress",
                      na.rm = TRUE
            ),
            HCo = sum(obs_by_reef$dist_type == "Heat Stress, CoTS Outbreak",
                      na.rm = TRUE
            ),
            Other = 0
        )
        event_counts_rf[["Other"]] <- event_counts_rf[["Other"]] +
            sum(!is.na(obs_by_reef$dist_type)) -
            sum(event_counts_rf)
        event_counts <- event_counts + event_counts_rf
        reef_unknown <- sum(grepl("Unknown", obs_by_reef$recov_year), na.rm = TRUE)
        overall_unknown_count <- overall_unknown_count + reef_unknown
        
        ## Get latitude & longitude
        coordinates <- st_coordinates(obs_by_reef$geometry)
        latitude <- coordinates[1, 2]
        longitude <- coordinates[1, 1]
        
        ## Get number of total disturbances
        num_total <- sum(obs_by_reef$is_disturbed, na.rm = TRUE)
        
        ## Convert baseline_vals to string
        baseline_vals <- paste(baseline_vals, collapse = ", ")
        
        ## Add all relevant information to data frame
        new_row <- c(
            reef_name, # reef or reef & site name
            latitude, # latitude of reef/site
            longitude, # longitude of reef/site
            obs_by_reef$AREA_DESCR[1], # management region
            num_total, # total number of disturbances in obs period
            num_single, # no. of singular dists in obs period
            prob_s_dist, # probability of single disturbance per year
            prob_s_impact, # probability of impact given single disturbance
            prob_s_recov, # recovery rate following single dist impact
            num_comp, # no. of cumulative dists in obs period
            prob_c_dist, # probability of compound disturbance per year
            prob_c_impact, # probability of impact given compound disturbance
            prob_c_recov, # recovery rate following compound dist impact
            reef_unknown, # number of unknown recovery times
            baseline_inferred, # if baseline inferred
            baseline_vals # baselines used in reef
        ) # number of unknown recovery times
        reef_df[nrow(reef_df) + 1, ] <- new_row
        # Replace all_reefs_sf with new info
        all_reefs_sf[all_reefs_sf$REEF_NAME == reef_name, ] <- obs_by_reef
    }
    
    # Return a list of the event counts and reef information dataframes
    return(list(event_counts = event_counts, all_reefs_sf = all_reefs_sf, reef_df = reef_df, overall_unknown_count = overall_unknown_count))
}