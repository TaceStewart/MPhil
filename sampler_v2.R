# Sample v2 for taking samples from a large df
# Sampler function commented
samplerv2 <- function(reef_df, sector_boundaries, sample_reefs, num_samples,
                      is_time_based, recov_th, recov_yrs, cots_dist, cyc_dist,
                      dhw_dist) {
  # Define the number of reefs per management area
  unique_mgmts <- unique(reef_df$sector)
  reefs_per_mgmt <- floor(num_samples / length(unique_mgmts))
  num_samples <- reefs_per_mgmt * length(unique_mgmts)
  
  # Take reefs_per_mgmt rows randomly from sample_reefs for each management area
  sample_reefs_all_yrs <- data.frame(matrix(ncol = ncol(sample_reefs),
                                            nrow = 0
  )) %>%
    setNames(names(sample_reefs))
  for (mgmt in unique_mgmts) {
    # Get the unique point ids for the management area
    mgmt_point_ids <- unique(sample_reefs[sample_reefs$sector == mgmt, "point_id"])
    
    # Sample reefs_per_mgmt reefs from mgmt_point_ids
    sample_mgmt_pts <- sample(mgmt_point_ids, reefs_per_mgmt, replace = FALSE)
    
    # Add the sampled reefs to the df
    sample_reefs_all_yrs <- rbind(
      sample_reefs_all_yrs,
      sample_reefs[sample_reefs$point_id %in% sample_mgmt_pts, ]
    )
  }
  sample_reefs_all_yrs$single_or_compound <- NA
  
  # Add a column for the disturbance type
  dist_names <- c("Wind Stress", "Heat Stress", "CoTS Outbreak")
  sample_reefs_all_yrs <- sample_reefs_all_yrs %>%
    mutate(dist_type = paste0(dist_names[Hs4MW_value >= cyc_dist &
                                           annMaxDHW_value >= dhw_dist &
                                           COTS_value >= cots_dist], collapse = ", "))
  
  # Make a df to store all sample reef info in
  sample_reefs_df <- unique(
    sample_reefs_all_yrs[, c("point_id", "sector", "point_loc")]
  ) %>%
    mutate(
      num_dist = NA,
      num_s_dist = NA,
      num_c_dist = NA,
      d_tot = NA,
      prob_s_dist = NA,
      prob_c_dist = NA,
      prob_s_impact = NA,
      prob_c_impact = NA,
      prob_s_recov = NA,
      prob_c_recov = NA
    )
  
  # Sample probabilities and recovery rates from those observed in each mgmt area
  for (row in 1:length(unique_mgmts)) {
    mgmt_index <- which(sector_boundaries$AREA_DESCR == unique_mgmts[row])
    sample_mgmt_indx <- which(sample_reefs_df$sector == unique_mgmts[row])
    
    # Sample the probability of impact from single disturbance for the management area
    reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    if (!any(!is.na(reef_df$prob_s_impact[reef_indx]))) {
      # If there are no single impact rates for the management area, use the
      # impact rates from the neighbouring management areas
      if (row == 1) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row + 1])
      } else if (row == length(unique_mgmts)) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1])
      } else {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1:row + 1])
      }
      reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    }
    mgmt_prob <- rnorm(n = reefs_per_mgmt,
                          mean = mean(reef_df$prob_s_impact[reef_indx],
                                      na.rm = TRUE),
                          sd = sd(reef_df$prob_s_impact[reef_indx],
                                  na.rm = TRUE))
    mgmt_prob[mgmt_prob < 0] <- 0
    mgmt_prob[mgmt_prob > 1] <- 1
    sample_reefs_df$prob_s_impact[sample_mgmt_indx] <- mgmt_prob
    
    # Sample the probability of impact from compound disturbance for the management area
    mgmt_index <- which(sector_boundaries$AREA_DESCR == unique_mgmts[row])
    reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    if (!any(!is.na(reef_df$prob_c_impact[reef_indx]))) {
      # If there are no compound impact rates for the management area, use the
      # impact rates from the neighbouring management areas
      if (row == 1) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row + 1])
      } else if (row == length(unique_mgmts)) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1])
      } else {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1:row + 1])
      }
      reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    }
    mgmt_prob <- rnorm(n = reefs_per_mgmt,
                          mean = mean(reef_df$prob_c_impact[reef_indx],
                                      na.rm = TRUE),
                          sd = sd(reef_df$prob_c_impact[reef_indx],
                                  na.rm = TRUE))
    mgmt_prob[mgmt_prob < 0] <- 0
    mgmt_prob[mgmt_prob > 1] <- 1
    sample_reefs_df$prob_c_impact[sample_mgmt_indx] <- mgmt_prob
    
    # Sample the probability of recovery from single disturbance for the management area
    mgmt_index <- which(sector_boundaries$AREA_DESCR == unique_mgmts[row])
    reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    if (!any(!is.na(reef_df$prob_s_recov[reef_indx]))) {
      # If there are no single recovery rates for the management area, use the
      # recovery rates from the neighbouring management areas
      if (row == 1) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row + 1])
      } else if (row == length(unique_mgmts)) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1])
      } else {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1:row + 1])
      }
      reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    }
    mgmt_prob <- rnorm(n = reefs_per_mgmt,
                          mean = mean(reef_df$prob_s_recov[reef_indx],
                                      na.rm = TRUE),
                          sd = sd(reef_df$prob_s_recov[reef_indx],
                                  na.rm = TRUE))
    mgmt_prob[mgmt_prob < 0] <- 0
    mgmt_prob[mgmt_prob > 1] <- 1
    sample_reefs_df$prob_s_recov[sample_mgmt_indx] <- mgmt_prob
    
    # Sample the probability of recovery from compound disturbance for the management area
    mgmt_index <- which(sector_boundaries$AREA_DESCR == unique_mgmts[row])
    reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    if (!any(!is.na(reef_df$prob_c_recov[reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index]]))) {
      # If there are no compound recovery rates for the management area, use the
      # recovery rates from the neighbouring management areas
      if (row == 1) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row + 1])
      } else if (row == length(unique_mgmts)) {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1])
      } else {
        mgmt_index <- which(sector_boundaries$AREA_DESCR %in% unique_mgmts[row - 1:row + 1])
      }
      reef_indx <- which(reef_df$sector %in% sector_boundaries$AREA_DESCR[mgmt_index])
    }
    mgmt_prob <- rnorm(n = reefs_per_mgmt,
                          mean = mean(reef_df$prob_c_recov[reef_indx],
                                      na.rm = TRUE),
                          sd = sd(reef_df$prob_c_recov[reef_indx],
                                  na.rm = TRUE))
    mgmt_prob[mgmt_prob < 0] <- 0
    mgmt_prob[mgmt_prob > 1] <- 1
    sample_reefs_df$prob_c_recov[sample_mgmt_indx] <- mgmt_prob
  }
  
  # For each unique sample reef
  for (reef in seq_len(length(sample_reefs_df$point_id))) {
    # Get the observations for that reef
    reef_obs <- sample_reefs_all_yrs[sample_reefs_all_yrs$point_id == sample_reefs_df$point_id[reef], ]
    
    # Calculate the probability of disturbance for the reef
    sample_reefs_df$d_tot[reef] <- sum(reef_obs$is_disturbed,
                                       na.rm = TRUE
    ) / nrow(reef_obs)
    
    # Determine if each disturbance at the reef is single or compound
    skip_next <- 0
    for (dist in seq_len(length(reef_obs$is_disturbed))) {
      if (skip_next > 0) {
        skip_next <- skip_next - 1
        next
      }
      
      # If there is a disturbance that year
      if (reef_obs$is_disturbed[dist]) {
        # If it is the final disturbance, it must be single
        if (dist == nrow(reef_obs)) {
          reef_obs$single_or_compound[dist] <- "Single"
          next
        }
        
        # Is there another disturbance at the reef within 1/r_single years?
        if (sum(
          reef_obs$is_disturbed[(dist + 1):nrow(reef_obs)] &
          reef_obs$year[(dist + 1):nrow(reef_obs)] <=
          reef_obs$year[dist] + 1 / as.numeric(sample_reefs_df$prob_s_recov[reef]),
          na.rm = TRUE
        ) > 0) {
          reef_obs$single_or_compound[dist] <- "Compound"
          skip_next <- sum(
            reef_obs$is_disturbed[(dist + 1):nrow(reef_obs)] &
              reef_obs$year[(dist + 1):nrow(reef_obs)] <=
              reef_obs$year[dist] + 1 / as.numeric(sample_reefs_df$prob_c_recov[reef]),
            na.rm = TRUE
          )
          
          # Combine the disturbance types for the cumulative disturbances
          reef_obs$dist_type[dist] <- paste0(
            reef_obs$dist_type[dist:(dist + skip_next)],
            collapse = ", "
          )
          
          # Remove future disturbance types from the cumulative disturbances
          reef_obs$dist_type[(dist + 1):(dist + skip_next)] <- NA
        } else {
          reef_obs$single_or_compound[dist] <- "Single"
        }
      }
      
      # Calculate the probability of single disturbance for the reef
      sample_reefs_df$prob_s_dist[reef] <- sum(
        reef_obs$is_disturbed &
          reef_obs$single_or_compound == "Single",
        na.rm = TRUE
      ) / nrow(reef_obs)
      
      # Calculate the probability of compound disturbance for the reef
      sample_reefs_df$prob_c_dist[reef] <- sum(
        reef_obs$is_disturbed &
          reef_obs$single_or_compound == "Compound",
        na.rm = TRUE
      ) / nrow(reef_obs)
      
      # Calculate number of disturbances at the reef
      sample_reefs_df$num_dist[reef] <- sum(
        reef_obs$is_disturbed,
        na.rm = TRUE
      )
      
      # Calculate the number of single disturbances at the reef
      sample_reefs_df$num_s_dist[reef] <- sum(
        reef_obs$is_disturbed &
          reef_obs$single_or_compound == "Single",
        na.rm = TRUE
      )
      
      # Calculate the number of cumulative disturbances at the reef
      sample_reefs_df$num_c_dist[reef] <- sum(
        reef_obs$is_disturbed &
          reef_obs$single_or_compound == "Compound",
        na.rm = TRUE
      )
    }
  }
  
  return(sample_reefs_df)
}
