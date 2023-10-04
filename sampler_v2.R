# Sample v2 for taking samples from a large df
# Sampler function commented
sampler <- function(reef_df, sector_boundaries, sample_reefs, num_samples, 
                    is_time_based, recov_th, recov_yrs, cots_dist, cyc_dist, 
                    dhw_dist) {
  # Define the number of reefs per management area
  unique_mgmts <- unique(reef_df$management_region)
  reefs_per_mgmt <- floor(num_samples / length(unique_mgmts))
  num_samples <- reefs_per_mgmt * length(unique_mgmts)
  
  # Take num_samples rows randomly from sample_reefs
  sample_rf_pt_ids <- sample(unique(sample_reefs$point_id), num_samples)
  sample_reefs_all_yrs <- sample_reefs[which(sample_reefs$point_id %in% sample_rf_pt_ids), ]
  
  # Make a df to store all sample reef info in
  sample_reefs_df <- data.frame(
    point_id = unique(sample_reefs_all_yrs$point_id),
    d_single = 0,
    d_comp = 0,
    r_single = 0,
    r_comp = 0
  )
  
  # Sample recovery rates from those observed in each mgmt area
  for (row in 1:length(unique_mgmts)) {
    mgmt_index <- sector_boundaries$AREA_DESCR == unique_mgmts[row]
    from_row <- (row - 1) * reefs_per_mgmt + 1
    to_row <- row * reefs_per_mgmt
    sample_reefs_df$r_single[from_row:to_row] <- sample(
      reef_df$r_single[!is.na(reef_df$r_single) &
                         reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]],
      reefs_per_mgmt,
      replace = TRUE
    )
    
    if (any(!is.na(reef_df$r_comp[reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]]))) {
      sample_reefs_df$r_comp[from_row:to_row] <- sample(
        reef_df$r_comp[!is.na(reef_df$r_comp) &
                         reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]],
        reefs_per_mgmt,
        replace = TRUE
      )
    } else { # Sampler needs to deal with unknowns better. Also, should sampling be from the data in the first place?
      mgmt_index <- sector_boundaries$AREA_DESCR == unique_mgmts[row - 1]
      from_row <- (row - 2) * reefs_per_mgmt + 1
      to_row <- (row - 1) * reefs_per_mgmt
      
      sample_reefs_df$r_comp[from_row:to_row] <- sample(
        reef_df$r_comp[!is.na(reef_df$r_comp) &
                         reef_df$sector == sector_boundaries$AREA_DESCR[mgmt_index]],
        reefs_per_mgmt,
        replace = TRUE
      )
    }
  }
  
  
  
  return(sample_reefs_df)
}
