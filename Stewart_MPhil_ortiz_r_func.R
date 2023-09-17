# For a given reef, this function finds the probability of single and compound 
# cyclones, and the probability of coral recovery following these events.
# Definition of compound: A disturbance occurs before coral has recovered from 
# the previous disturbance
ortiz_r_func <- function(obs_by_reef, cots_dist, cyc_dist, dhw_dist) {
  # Set epsilon
  eps <- 0.05 # 5% 
  
  # If there are no disturbances
  if (max(obs_by_reef$COTS_value) < cots_dist & 
      max(obs_by_reef$Hs4MW_value) < cyc_dist &
      max(obs_by_reef$DHW_value) < dhw_dist) {
    
    
  } else {
    # Loop through each disturbance
    distIndices <- which(!is.na(obs_by_reef$recovYear))
    for (distIndex in distIndices) {
      # If the recovery time is unknown
      if (!is.na(obs_by_reef$recovYear[distIndex]) &
          grepl("Unknown", obs_by_reef$recovYear[distIndex])) { 
        obs_by_reef$ortizRecov[distIndex] <- "Unknown - no recovery time"
      } else {
        # Calculate Ortiz growth rate
        priorCover <- obs_by_reef$COVER[distIndex]
        postRecovCover <- obs_by_reef$COVER[obs_by_reef$YEAR == obs_by_reef$recovYear[distIndex]][1]
        
        # If the cover in year of recovery exists
        if (length(postRecovCover) > 0 ) {
          N_T <- (1 + eps)*postRecovCover                # Final coral value
          N_t <- (1 + eps)*priorCover                    # Initial coral value
          T_minus_t <- obs_by_reef$recovTime[distIndex]  # Time to recover
          obs_by_reef$ortizRecov[distIndex] <- (log(N_T) - log(N_t))/T_minus_t
          
        } else { # If the cover in year of recovery does not exist
          # Check if there's any sign of growth in post observations
          smallestPostObsIndex <- which.min(obs_by_reef$COVER[(distIndex + 1):nrow(obs_by_reef)]) + distIndex
          if (!is.na(smallestPostObsIndex) &
              smallestPostObsIndex != nrow(obs_by_reef)) {
            # Get the next obs greater than that
            growthObsIndex <- which(obs_by_reef$COVER[(smallestPostObsIndex + 1):nrow(obs_by_reef)] 
                                    > obs_by_reef$COVER[smallestPostObsIndex])[1] + smallestPostObsIndex
            
            N_T <- (1 + eps)*obs_by_reef$COVER[growthObsIndex]        # Final coral value
            N_t <- (1 + eps)*obs_by_reef$COVER[smallestPostObsIndex]  # Initial coral value
            T_minus_t <- obs_by_reef$YEAR[growthObsIndex] - 
              obs_by_reef$YEAR[smallestPostObsIndex]                  # Time to recover
            obs_by_reef$ortizRecov[distIndex] <- (log(N_T) - log(N_t))/T_minus_t
            
          } else { # If there's no sign of growth in post observations
            obs_by_reef$ortizRecov[distIndex] <- "Unknown - no growth post disturbance"
          }
        }
        
      }
    }
  }
  
  # Return modified obs_by_reef df
  return(obs_by_reef)
}