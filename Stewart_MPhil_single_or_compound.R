# For a given reef, this function finds the probability of single and compound
# cyclones, and the probability of coral recovery following these events.
# Definition of compound if is_time_based is FALSE: i.e., recovery based,
# A disturbance occurs before coral has recovered from the previous disturbance
# Definition of compound if is_time_based is TRUE: A disturbance occurs within
# <recov_yrs> of a previous disturbance
single_or_compound <- function(obs_by_reef, is_time_based, recov_th, recov_yrs,
                               cots_dist, cyc_dist, dhw_dist,
                               infer_baseline = 0, epsilon = 0.05,
                               baseline_str = "mean") {
  # Create vector of disturbance string names
  dist_names <- c("CoTS Outbreak", "Wind Stress", "Heat Stress")
  baseline_inferred <- NA
  baseline_vals <- c()
  
  # If there are no disturbances
  if (max(obs_by_reef$COTS_value) < cots_dist &&
      max(obs_by_reef$Hs4MW_value) < cyc_dist &&
      max(obs_by_reef$DHW_value) < dhw_dist) {
  } else { # If there are disturbances
    # Get indices of disturbance years
    dist_indices <- which(obs_by_reef$is_disturbed)
    
    # If time based
    if (is_time_based) {
      for (dist_index in dist_indices) {
        # If there are any disturbances between this one and recov_yrs from now
        current_yr <- obs_by_reef$YEAR[dist_index]
        if (any(obs_by_reef$is_disturbed[obs_by_reef$YEAR > current_yr &
                                         obs_by_reef$YEAR < current_yr + recov_yrs])) {
          # Get a list of those disturbance index/es
          next_dist_s <- which(obs_by_reef$is_disturbed[obs_by_reef$YEAR > current_yr &
                                                          obs_by_reef$YEAR < current_yr + recov_yrs]) + dist_index
          skip_next <- length(next_dist_s)
          
          # While there are still disturbances within recov_yrs of the final dist in compound cluster
          final_dist_index <- next_dist_s[length(next_dist_s)]
          current_yr <- obs_by_reef$YEAR[final_dist_index]
          while (any(obs_by_reef$is_disturbed[obs_by_reef$YEAR > current_yr &
                                              obs_by_reef$YEAR < current_yr + recov_yrs])) {
            # Get a list of those disturbance index/es
            next_dist_s <- which(obs_by_reef$is_disturbed[obs_by_reef$YEAR > current_yr &
                                                            obs_by_reef$YEAR < current_yr + recov_yrs]) + final_dist_index
            skip_next <- skip_next + length(next_dist_s)
            final_dist_index <- next_dist_s[length(next_dist_s)]
            current_yr <- obs_by_reef$YEAR[final_dist_index]
          }
          
          # Recovery year is the next disturbance year + recovery time
          obs_by_reef$recov_year[dist_index] <- max(obs_by_reef$YEAR[next_dist_s]) + 
            recov_yrs
        } else {
          # Assume recovery year is the disturbance year + recovery time
          obs_by_reef$recov_year[dist_index] <- obs_by_reef$YEAR[dist_index] + recov_yrs
        }
      }
    } else { # If recovery based
      # Set first baseline as the max coral value before first disturbance
      max_pre_dist_cover <- ifelse(dist_indices[1] == 1,
                                   NA,
                                   max(obs_by_reef$COVER[1:dist_indices[1] - 1])
      )
      
      # If coral cover at first disturbance is below max_pre_dist_cover*(1 - epsilon)
      if (!is.na(max_pre_dist_cover) &&
          obs_by_reef$COVER[dist_indices[1]] < max_pre_dist_cover * (1 - epsilon)) {
        # Then we have a reef with a baseline
        baseline <- max_pre_dist_cover * (1 - epsilon)
        baseline_vals <- c(baseline_vals, baseline)
        baseline_inferred <- FALSE
      } else if (infer_baseline) {
        # Find local maxima within the epsilon range
        local_maxima <- findpeaks(obs_by_reef$COVER, nups = 1, ndowns = 1)
        maxima_values <- local_maxima[, 1]
        
        if (baseline_str == "min") {
          # Calculate the minimum of the local maxima
          baseline <- ifelse(length(local_maxima) > 0,
                             min(maxima_values),
                             min(obs_by_reef$COVER)
          )
        } else if (baseline_str == "mean") {
          # Calculate the average of the local maxima
          baseline <- ifelse(length(local_maxima) > 0,
                             mean(maxima_values),
                             mean(obs_by_reef$COVER)
          )
        } else if (baseline_str == "max") {
          # Calculate the max of the local maxima
          baseline <- ifelse(length(local_maxima) > 0,
                             max(maxima_values),
                             max(obs_by_reef$COVER)
          )
        }
        
        baseline_vals <- c(baseline_vals, baseline)
        baseline_inferred <- TRUE
      } else {
        # Then this reef does not have a baseline
        baseline <- NA
        baseline_inferred <- FALSE
      }
      # Initialise skip variable
      skip_next <- 0
      if (!is.na(baseline)) {
        # For each disturbance
        # dist_index <- dist_indices[1] # for testing
        for (dist_index in dist_indices) {
          if (dist_index > 2 &&
              is.na(max_pre_dist_cover)) {
            max_pre_dist_cover <- max(obs_by_reef$COVER[1:dist_index - 1])
          }
          
          # Check if there's a new baseline
          if (dist_index > 2 &&
              max(obs_by_reef$COVER[1:dist_index - 1]) > max_pre_dist_cover) {
            max_pre_dist_cover <- max(obs_by_reef$COVER[1:dist_index - 1])
            baseline <- max_pre_dist_cover * (1 - epsilon)
            baseline_vals <- c(baseline_vals, baseline)
          }
          # If we need to skip this disturbance as it was already counted in the last one
          if (skip_next > 0) {
            # Take one from the skip variable and move on
            skip_next <- skip_next - 1
          } else {
            # Add disturbance type
            obs_by_reef$dist_type[dist_index] <- paste0(
              dist_names[c(
                obs_by_reef$COTS_value[dist_index] >= cots_dist,
                obs_by_reef$Hs4MW_value[dist_index] >= cyc_dist,
                obs_by_reef$DHW_value[dist_index] >= dhw_dist
              )],
              collapse = ", "
            )
            # Determine if reef is impacted & log in df
            obs_by_reef$is_impacted[dist_index] <- obs_by_reef$COVER[dist_index] < baseline
            
            # If it is impacted, find recovery year
            if (obs_by_reef$is_impacted[dist_index]) {
              # If we're in the last row, we have no post obs
              if (dist_index == nrow(obs_by_reef)) {
                # Set recovery year to unknown
                obs_by_reef$recov_year[dist_index] <- "Unknown - no post obs"
              } else {
                # Find next observation with at least max_pre_dist_cover * (1 - recov_th)
                post_rec_index <- which(obs_by_reef$COVER[dist_index + 1:nrow(obs_by_reef)] > 
                                          max_pre_dist_cover * (1 - recov_th) & 
                                          obs_by_reef$COVER[dist_index + 1:nrow(obs_by_reef)] >
                                          min(obs_by_reef$COVER[dist_index:nrow(obs_by_reef) - 1]))[1] + dist_index
                # If there are no obs with at least (1 - recov_th) * max_pre_dist_cover
                if (is.na(post_rec_index)) {
                  obs_by_reef$recov_year[dist_index] <- "Unknown - no post-disturbance growth"
                  
                  # If there are any disturbances afterwards
                  if (any(obs_by_reef$is_disturbed[(dist_index + 1):nrow(obs_by_reef)])) {
                    # Pair these disturbances with the first instance and make it compound
                    next_dists <- which(obs_by_reef$is_disturbed[(dist_index + 1):nrow(obs_by_reef)]) + dist_index
                    dist_type <- obs_by_reef$dist_type[dist_index]
                    for (event in next_dists) {
                      event_disturb <- dist_names[c(
                        obs_by_reef$COTS_value[event] >= cots_dist,
                        obs_by_reef$Hs4MW_value[event] >= cyc_dist,
                        obs_by_reef$DHW_value[event] >= dhw_dist
                      )]
                      if (length(event_disturb > 1)) {
                        event_disturb <- paste(event_disturb, collapse = ", ")
                      }
                      dist_type <- paste0(
                        c(
                          dist_type,
                          event_disturb
                        ),
                        collapse = ", "
                      )
                      obs_by_reef$dist_type[event] <- NA
                      obs_by_reef$recov_year[event] <- NA
                    }
                    obs_by_reef$dist_type[dist_index] <- dist_type
                    
                    # Make sure those dists are skipped
                    skip_next <- length(next_dists)
                  }
                } else { # If there are obs with at least baseline
                  # If there are any disturbances between disturbance and recovery
                  if (dist_index + 1 <= post_rec_index - 1 &&
                      any(obs_by_reef$is_disturbed[(dist_index + 1):(post_rec_index - 1)])) {
                    # Pair these disturbances with the first instance and make it compound
                    next_dists <- which(obs_by_reef$is_disturbed[(dist_index + 1):(post_rec_index - 1)]) + dist_index
                    dist_type <- obs_by_reef$dist_type[dist_index]
                    for (event in next_dists) {
                      event_disturb <- dist_names[c(
                        obs_by_reef$COTS_value[event] >= cots_dist,
                        obs_by_reef$Hs4MW_value[event] >= cyc_dist,
                        obs_by_reef$DHW_value[event] >= dhw_dist
                      )]
                      if (length(event_disturb > 1)) {
                        event_disturb <- paste(event_disturb, collapse = ", ")
                      }
                      dist_type <- paste0(
                        c(
                          dist_type,
                          event_disturb
                        ),
                        collapse = ", "
                      )
                      obs_by_reef$dist_type[event] <- NA
                      obs_by_reef$recov_year[event] <- NA
                    }
                    obs_by_reef$dist_type[dist_index] <- dist_type
                    
                    # Make sure those dists are skipped
                    skip_next <- length(next_dists)
                  }
                  obs_by_reef$recov_year[dist_index] <- obs_by_reef$YEAR[post_rec_index]
                }
              }
            }
          }
        }
      }
    }
  }
  # List of indices of disturbances with recovYr so we don't loop through all rows
  dist_indices <- which(!is.na(obs_by_reef$recov_year) &
                          (!grepl("Unknown", obs_by_reef$recov_year)))
  
  # If there are disturbances with recovYr
  if (length(dist_indices) > 0) {
    for (index in dist_indices) { # index <- dist_indices[5] # for testing
      if (any(
        obs_by_reef$YEAR > obs_by_reef$YEAR[index] & #  If any coral obs until recovery year exist
        obs_by_reef$YEAR < obs_by_reef$recov_year[index],
        na.rm = TRUE
      ) &&
      any(
        obs_by_reef$is_disturbed[obs_by_reef$YEAR > obs_by_reef$YEAR[index] & # And if there is another disturbance between this one and recovery
                                 obs_by_reef$YEAR < obs_by_reef$recov_year[index]],
        na.rm = TRUE
      )) {
        # Get index of all following disturbances until recovery
        nextdist_indices <- which(obs_by_reef$is_disturbed[obs_by_reef$YEAR > obs_by_reef$YEAR[index] &
                                                             obs_by_reef$YEAR < obs_by_reef$recov_year[index]]) + index
        
        # Add disturbance type
        dist_type <- obs_by_reef$dist_type[index]
        for (event in nextdist_indices) {
          event_disturb <- dist_names[c(
            obs_by_reef$COTS_value[event] >= cots_dist,
            obs_by_reef$Hs4MW_value[event] >= cyc_dist,
            obs_by_reef$DHW_value[event] >= dhw_dist
          )]
          if (length(event_disturb > 1)) {
            event_disturb <- paste(event_disturb, collapse = ", ")
          }
          dist_type <- paste0(
            c(
              dist_type,
              event_disturb
            ),
            collapse = ", "
          )
          obs_by_reef$dist_type[event] <- NA
          obs_by_reef$recov_year[event] <- NA
        }
        obs_by_reef$dist_type[index] <- dist_type
      } else {
        # Get disturbance/s for that year (there may be multiple types in one year)
        dist_type <- dist_names[c(
          obs_by_reef$COTS_value[index] >= cots_dist,
          obs_by_reef$Hs4MW_value[index] >= cyc_dist,
          obs_by_reef$DHW_value[index] >= dhw_dist
        )]
        if (length(dist_type > 1)) {
          dist_type <- paste(dist_type, collapse = ", ")
        }
        
        # Add disturbance type
        obs_by_reef$dist_type[index] <- dist_type
      } # end of single/compound if-else statement
      
      # Add recovery time
      obs_by_reef$recov_time[index] <- as.numeric(obs_by_reef$recov_year[index]) -
        obs_by_reef$YEAR[index]
      if (is_time_based) {
        obs_by_reef$is_impacted[index] <- runif(1)
        obs_by_reef$r_given_impact[index] <- 1 / recov_yrs
      } else {
        obs_by_reef$r_given_impact[index] <- ifelse(obs_by_reef$is_impacted[index],
                                                    1 / (obs_by_reef$recov_time[index]),
                                                    NA
        )
        
      }
      
      # If Inf was produced, there was a 0 year recovery time, so r is 100%
      if (obs_by_reef$r_given_impact[index] == Inf) {
        obs_by_reef$r_given_impact[index] <- 1
      }
    } # end of disturbance for loop
  }
  
  for (disturbance in which(!is.na(obs_by_reef$dist_type))) {
    # If there are any commas, it's compound
    if (grepl(",", obs_by_reef$dist_type[disturbance])) {
      obs_by_reef$single_or_compound[disturbance] <- "Compound"
    } else {
      # If not, it's single
      obs_by_reef$single_or_compound[disturbance] <- "Single"
    }
  }
  
  # Return recovery and disturbance probabilities
  return(list(obs_by_reef, baseline_inferred, baseline_vals)) 
}
