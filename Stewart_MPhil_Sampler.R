sampler <- function(reef_df, sector_boundaries, num_samples, isTimeBased, 
                    recov_th, recov_yrs, cots_dist, cyc_dist, dhw_dist) {
  # Define the number of reefs per management area
  unique_mgmts <- unique(reef_df$sector)
  reefs_per_mgmt <- floor(num_samples / length(unique_mgmts))
  num_samples <- reefs_per_mgmt * length(unique_mgmts)
  
  # Initialise a list of random points
  random_points <- vector("integer", length = num_samples)
  
  # Get random locations in each management area for sample reefs 
  for (row in 1:length(unique_mgmts)) {
    mgmt_index <- sector_boundaries$AREA_DESCR == unique_mgmts[row]
    from_row <- (row - 1)*reefs_per_mgmt + 1
    to_row <- row*reefs_per_mgmt
    random_points[from_row:to_row] <- sector_boundaries$geometry[mgmt_index] %>%
      st_sample(size = reefs_per_mgmt, 
                type = "random")
  }
  
  # Extract coordinates and create a data frame
  coordinates <- lapply(random_points, sf::st_coordinates)
  sample_reefs_df <- data.frame(x = sapply(coordinates, function(x) x[, 1]),
                                y = sapply(coordinates, function(x) x[, 2]))
  
  # Rename columns
  colnames(sample_reefs_df) <- c("LONGITUDE", "LATITUDE")
  
  # Assign hexagon IDs and sample values
  sample_reefs_df <- sample_reefs_df %>%
    mutate(point_id = row_number(),
           d_single = 0,
           d_comp = 0,
           r_single = 0,
           r_comp = 0)
  
  # Sample probabilities and recovery rates from those observed in each mgmt area
  for (row in 1:length(unique_mgmts)) {
    mgmt_index <- sector_boundaries$AREA_DESCR == unique_mgmts[row]
    from_row <- (row - 1)*reefs_per_mgmt + 1
    to_row <- row*reefs_per_mgmt
    sample_reefs_df$r_single[from_row:to_row] <- sample(reef_df$r_single[!is.na(reef_df$r_single) &
                                                                           reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]],
                                                        reefs_per_mgmt,
                                                        replace = TRUE)    
    
    if (any(!is.na(reef_df$r_comp[reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]]))) {
      sample_reefs_df$r_comp[from_row:to_row] <- sample(reef_df$r_comp[!is.na(reef_df$r_comp) &
                                                                         reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]],
                                                        reefs_per_mgmt,
                                                        replace = TRUE)
    } else { # Sampler needs to deal with unknowns better. Also, should sampling be from the data in the first place?
      mgmt_index <- sector_boundaries$AREA_DESCR == unique_mgmts[row-1]
      from_row <- (row - 2)*reefs_per_mgmt + 1
      to_row <- (row-1)*reefs_per_mgmt
      
      sample_reefs_df$r_comp[from_row:to_row] <- sample(reef_df$r_comp[!is.na(reef_df$r_comp) &
                                                                         reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]],
                                                        reefs_per_mgmt,
                                                        replace = TRUE)
    }
  }
  
  # Find closest disturbance observations
  sample_reefs_df <- dist_finder(sample_reefs_df, isTimeBased, recov_th, 
                                 recov_yrs, cots_dist, cyc_dist, dhw_dist)
  
  return(sample_reefs_df)
}
