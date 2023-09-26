# For a given reef, this function finds the probability of single and compound 
# cyclones, and the probability of coral recovery following these events.
# Definition of compound if isTimeBased is FALSE: i.e., recovery based, 
# A disturbance occurs before coral has recovered from the previous disturbance
# Definition of compound if isTimeBased is TRUE: A disturbance occurs within 
# <recov_yrs> of a previous disturbance
single_or_compound <- function(obs_by_reef, isTimeBased, recov_th, recov_yrs,
                               cots_dist, cyc_dist, dhw_dist, inferBaseline = 0,
                               epsilon = 0.05, baseline_str = "mean") {
  
  # Create vector of disturbance string names
  distNames <- c("CoTS Outbreak", "Wind Stress", "Heat Stress")
  
  # If there are no disturbances
  if (max(obs_by_reef$COTS_value) < cots_dist & 
      max(obs_by_reef$Hs4MW_value) < cyc_dist &
      max(obs_by_reef$DHW_value) < dhw_dist) {
    
  } else { # If there are disturbances 
    # Get indices of disturbance years
    distIndices <- which(obs_by_reef$isDisturbed)
    
    if(distIndices[1] == 1) {
      
    } else {
      # Set first baseline as the max coral value before first disturbance
      maxPreDistCover <- max(obs_by_reef$COVER[1:distIndices[1] - 1]) 
      
      # If coral cover at first disturbance is below maxPreDistCover*(1 - epsilon) 
      if(obs_by_reef$COVER[distIndices[1]] < maxPreDistCover*(1 - epsilon)) {
        # Then we have a reef with a baseline
        baseline <- maxPreDistCover*(1 - epsilon)
      } else if (inferBaseline) { 
        # Find local maxima within the epsilon range
        local_maxima <- findpeaks(obs_by_reef$COVER, nups=1, ndowns=1)
        maxima_values <- local_maxima[,1]
        
        if (baseline_str == "min") {
          # Calculate the minimum of the local maxima
          baseline <- ifelse(length(local_maxima) > 0,
                             min(maxima_values),
                             min(obs_by_reef$COVER))
        } else if (baseline_str == "mean") {
          # Calculate the average of the local maxima
          baseline <- ifelse(length(local_maxima) > 0,
                             mean(maxima_values),
                             mean(obs_by_reef$COVER))
        } else if (baseline_str == "max") {
          # Calculate the max of the local maxima
          baseline <- ifelse(length(local_maxima) > 0,
                             max(maxima_values),
                             max(obs_by_reef$COVER))
        }
        
      } else {
        # Then this reef does not have a baseline
        baseline <- NA
      }
      # Initialise skip variable
      skipNext <- 0
      if (!is.na(baseline)) {
        # For each disturbance
        # distIndex <- distIndices[1] # for testing
        for (distIndex in distIndices) {
          # Check if there's a new baseline
          if(max(obs_by_reef$COVER[1:distIndex - 1]) > maxPreDistCover) {
            maxPreDistCover <- max(obs_by_reef$COVER[1:distIndex - 1])
            baseline <- maxPreDistCover*(1 - epsilon)
          }
          
          # If we need to skip this disturbance as it was already counted in the last one
          if (skipNext > 0) {
            # Take one from the skip variable and move on
            skipNext <- skipNext - 1 
          } else { 
            # Add disturbance type
            obs_by_reef$distType[distIndex] <- paste0(distNames[c(obs_by_reef$COTS_value[distIndex] >= cots_dist, 
                                                                  obs_by_reef$Hs4MW_value[distIndex] >= cyc_dist, 
                                                                  obs_by_reef$DHW_value[distIndex] >= dhw_dist)],
                                                      collapse = ', ')
            # If time based
            if (isTimeBased) { 
              # If there are any disturbances between this one and recov_yrs from now
              currentYr <- obs_by_reef$YEAR[distIndex]
              if (any(obs_by_reef$isDisturbed[obs_by_reef$YEAR > currentYr &
                                              obs_by_reef$YEAR < currentYr + recov_yrs])) {   
                # Get a list of those disturbance index/es
                nextDist_s <- which(obs_by_reef$isDisturbed[obs_by_reef$YEAR > currentYr &
                                                              obs_by_reef$YEAR < currentYr + recov_yrs]) + distIndex
                skipNext <- length(nextDist_s)
                
                # While there are still disturbances within recov_yrs of the final dist in compound cluster
                finalDistIndex <- nextDist_s[length(nextDist_s)]
                currentYr <- obs_by_reef$YEAR[finalDistIndex]
                while (any(obs_by_reef$isDisturbed[obs_by_reef$YEAR > currentYr &
                                                   obs_by_reef$YEAR < currentYr + recov_yrs])) {
                  # Get a list of those disturbance index/es
                  nextDist_s <- which(obs_by_reef$isDisturbed[obs_by_reef$YEAR > currentYr &
                                                                obs_by_reef$YEAR < currentYr + recov_yrs]) + finalDistIndex
                  skipNext <- skipNext + length(nextDist_s)
                  finalDistIndex <- nextDist_s[length(nextDist_s)]
                  currentYr <- obs_by_reef$YEAR[finalDistIndex]
                }
                
                # Recovery year is the next disturbance year + recovery time
                obs_by_reef$recovYear[distIndex] <- max(obs_by_reef$YEAR[nextDist_s]) + recov_yrs
                
              } else {
                # Assume recovery year is the disturbance year + recovery time
                obs_by_reef$recovYear[distIndex] <- obs_by_reef$YEAR[distIndex] + recov_yrs
              }
            } else { # If recovery based
              # Determine if reef is impacted & log in df
              obs_by_reef$isImpacted[distIndex] <- obs_by_reef$COVER[distIndex] < baseline
              
              # If it is impacted, find recovery year
              if (obs_by_reef$isImpacted[distIndex]) {
                # If we're in the last row, we have no post obs
                if (distIndex == nrow(obs_by_reef)) {
                  # Set recovery year to unknown
                  obs_by_reef$recovYear[distIndex] <- "Unknown - no post obs"
                } else { 
                  # Find next observation with at least maxPreDistCover*recov_th
                  postRecIndex <- which(obs_by_reef$COVER[distIndex:nrow(obs_by_reef)] > maxPreDistCover*recov_th)[1] + distIndex - 1
                  # If there are no obs with at least recov_th*maxPreDistCover
                  if (is.na(postRecIndex)) {
                    obs_by_reef$recovYear[distIndex] <- "Unknown - no post-disturbance growth"
                    
                    # If there are any disturbances afterwards
                    if (any(obs_by_reef$isDisturbed[(distIndex + 1):nrow(obs_by_reef)])) {
                      # Pair these disturbances with the first instance and make it compound
                      nextDists <- which(obs_by_reef$isDisturbed[(distIndex + 1):nrow(obs_by_reef)]) + distIndex
                      distType <- obs_by_reef$distType[distIndex]
                      for(event in nextDists) {
                        eventDisturb <- distNames[c(obs_by_reef$COTS_value[event] >= cots_dist, 
                                                    obs_by_reef$Hs4MW_value[event] >= cyc_dist, 
                                                    obs_by_reef$DHW_value[event] >= dhw_dist)]
                        if (length(eventDisturb > 1)) {
                          eventDisturb <- paste(eventDisturb, collapse = ", ")
                        } 
                        distType <- paste0(c(distType,
                                             eventDisturb),
                                           collapse = ", ")
                        obs_by_reef$distType[event] <- NA
                        obs_by_reef$recovYear[event] <- NA
                      }
                      obs_by_reef$distType[distIndex] <- distType
                      
                      # Make sure those dists are skipped
                      skipNext <- length(nextDists) 
                    }
                  } else { # If there are obs with at least baseline
                    # If there are any disturbances between disturbance and recovery
                    if (distIndex + 1 <= postRecIndex - 1 &
                        any(obs_by_reef$isDisturbed[(distIndex + 1):(postRecIndex - 1)])) {
                      # Pair these disturbances with the first instance and make it compound
                      nextDists <- which(obs_by_reef$isDisturbed[(distIndex + 1):(postRecIndex - 1)]) + distIndex
                      distType <- obs_by_reef$distType[distIndex]
                      for(event in nextDists) {
                        eventDisturb <- distNames[c(obs_by_reef$COTS_value[event] >= cots_dist, 
                                                    obs_by_reef$Hs4MW_value[event] >= cyc_dist, 
                                                    obs_by_reef$DHW_value[event] >= dhw_dist)]
                        if (length(eventDisturb > 1)) {
                          eventDisturb <- paste(eventDisturb, collapse = ", ")
                        } 
                        distType <- paste0(c(distType,
                                             eventDisturb),
                                           collapse = ", ")
                        obs_by_reef$distType[event] <- NA
                        obs_by_reef$recovYear[event] <- NA
                      }
                      obs_by_reef$distType[distIndex] <- distType
                      
                      # Make sure those dists are skipped
                      skipNext <- length(nextDists)
                    }
                    obs_by_reef$recovYear[distIndex] <- obs_by_reef$YEAR[postRecIndex]
                  } 
                } 
              }
            } 
          }
        } 
      }
    }
  }  
  # List of indices of disturbances with recovYr so we don't loop through all rows
  distIndices <- which(!is.na(obs_by_reef$recovYear) & 
                         (!grepl("Unknown", obs_by_reef$recovYear)))
  
  # If there are disturbances with recovYr
  if (length(distIndices) > 0) {
    for (index in distIndices) { # index <- distIndices[5] # for testing
      if (any(obs_by_reef$YEAR > obs_by_reef$YEAR[index] & #  If any coral obs until recovery year exist
              obs_by_reef$YEAR < obs_by_reef$recovYear[index], 
              na.rm = TRUE) & 
          any(obs_by_reef$isDisturbed[obs_by_reef$YEAR > obs_by_reef$YEAR[index] &  # And if there is another disturbance between this one and recovery
                                      obs_by_reef$YEAR < obs_by_reef$recovYear[index]], 
              na.rm = TRUE)) {
        
        # Get index of all following disturbances until recovery
        nextDistIndices <- which(obs_by_reef$isDisturbed[obs_by_reef$YEAR > obs_by_reef$YEAR[index] & 
                                                           obs_by_reef$YEAR < obs_by_reef$recovYear[index]]) + index
        
        # Add disturbance type
        distType <- obs_by_reef$distType[index]
        for(event in nextDistIndices) {
          eventDisturb <- distNames[c(obs_by_reef$COTS_value[event] >= cots_dist, 
                                      obs_by_reef$Hs4MW_value[event] >= cyc_dist, 
                                      obs_by_reef$DHW_value[event] >= dhw_dist)]
          if (length(eventDisturb > 1)) {
            eventDisturb <- paste(eventDisturb, collapse = ", ")
          } 
          distType <- paste0(c(distType,
                               eventDisturb),
                             collapse = ", ")
          obs_by_reef$distType[event] <- NA
          obs_by_reef$recovYear[event] <- NA
        }
        obs_by_reef$distType[index] <- distType
        
      } else { 
        # Get disturbance/s for that year (there may be multiple types in one year)
        distType <- distNames[c(obs_by_reef$COTS_value[index] >= cots_dist, 
                                obs_by_reef$Hs4MW_value[index] >= cyc_dist, 
                                obs_by_reef$DHW_value[index] >= dhw_dist)]
        if (length(distType > 1)) {
          distType <- paste(distType, collapse = ", ")
        } 
        
        # Add disturbance type
        obs_by_reef$distType[index] <- distType
      } # end of single/compound if-else statement
      
      # Add recovery time
      obs_by_reef$recovTime[index] <- as.numeric(obs_by_reef$recovYear[index]) - 
        obs_by_reef$YEAR[index]
      obs_by_reef$r_given_impact[index] <- ifelse(obs_by_reef$isImpacted[index],
                                                  1/(obs_by_reef$recovTime[index]),
                                                  NA)
      if(obs_by_reef$r_given_impact[index] == Inf) {
        obs_by_reef$r_given_impact[index] <- 1
      }
    } # end of disturbance for loop
  } else { # If there are no disturbances with recovYr
    
  } # end of disturbance existence if-else statement
  
  
  for(disturbance in which(!is.na(obs_by_reef$distType))) {
    # If there are any commas, it's compound
    if(grepl(',', obs_by_reef$distType[disturbance])) {
      obs_by_reef$single_or_compound[disturbance] <- "Compound"
    } else {
      # If not, it's single
      obs_by_reef$single_or_compound[disturbance] <- "Single"
    }
  }
  
  # Return recovery and disturbance probabilities
  return(obs_by_reef)
} 