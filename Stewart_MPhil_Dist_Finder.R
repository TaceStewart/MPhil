dist_finder <- function(df, isTimeBased, recov_th, recov_yrs,
                        cots_dist, cyc_dist, dhw_dist) {
  
  # df <- sample_reefs_df #test
  # Load the shapefile
  shapefile_path <- paste0("../Data/GBRMPA/Management_Areas_of_the_Great_Barrier",
                           "_Reef_Marine_Park/Management_Areas_of_the_Great_",
                           "Barrier_Reef_Marine_Park.shp")
  sector_boundaries <- st_read(shapefile_path, quiet = TRUE)
  
  # Import Data
  cots_data_raw <- read.csv("../Data/GBRdata/CoTS_data.csv")
  cyc_data_raw <- read.csv("../Data/GBRdata/Cyclones_data.csv")
  dhw_data_raw <- read.csv("../Data/GBRdata/DHW_data.csv")
  
  # Create spatial points for the reefs using the longitude and latitude coordinates
  df <- st_as_sf(df, 
                 coords = c("LONGITUDE", "LATITUDE"), 
                 crs = st_crs(sector_boundaries))
  
  # Perform a spatial join to determine which area each reef belongs to
  df <- st_join(df, sector_boundaries[, "AREA_DESCR"])
  
  # Transform columns to rows, summarise cots at each Reef ID
  cots_data <- cots_data_raw %>%
    pivot_longer(cols = starts_with("COTS"),
                 names_to = "year",
                 values_to = "COTS_value",
                 values_drop_na = FALSE) %>%
    mutate(year = stringr::str_extract(year, "\\d{4}")) %>%
    group_by(REEF_ID,year) %>%
    summarize(COTS_value = sum(COTS_value, na.rm = TRUE)) %>%
    suppressMessages()
  
  cyc_data <- cyc_data_raw %>%
    pivot_longer(cols = starts_with("Hs4MW"),
                 names_to = "year",
                 values_to = "Hs4MW_value",
                 values_drop_na = FALSE) %>%
    mutate(year = stringr::str_extract(year, "\\d{4}")) %>%
    group_by(REEF_ID,year) %>%
    summarize(Hs4MW_value = mean(Hs4MW_value, na.rm = TRUE)) %>%
    suppressMessages()
  
  dhw_data <- dhw_data_raw %>%
    pivot_longer(cols = starts_with("annMaxDHW"),
                 names_to = "year",
                 values_to = "annMaxDHW_value",
                 values_drop_na = FALSE) %>%
    mutate(year = stringr::str_extract(year, "\\d{4}")) %>%
    group_by(REEF_ID,year) %>%
    summarize(annMaxDHW_value = mean(annMaxDHW_value, na.rm = TRUE)) %>%
    suppressMessages()
  
  # Transform long and lat into polygons for cots_data
  cots_sf_unique <- cots_data_raw %>%
    st_as_sf(coords = c("LONG", "LAT"), crs = st_crs(sector_boundaries)) %>%
    group_by(REEF_ID) %>%
    summarise(geometry = st_combine(geometry))
  
  cyc_sf_unique <- cyc_data_raw %>%
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(sector_boundaries)) %>%
    group_by(REEF_ID) %>%
    summarise(geometry = st_combine(geometry))
  
  dhw_sf_unique <- dhw_data_raw %>%
    st_as_sf(coords = c("lon", "lat"), crs = st_crs(sector_boundaries)) %>%
    group_by(REEF_ID) %>%
    summarise(geometry = st_combine(geometry))
  
  # Join unique geometry with cots obs
  cots_data_sf <- left_join(cots_sf_unique, cots_data) %>%
    dplyr::select(REEF_ID, year, COTS_value, geometry) %>%
    suppressMessages()
  cots_data_sf$year <- as.numeric(cots_data_sf$year)
  
  cyc_data_sf <- left_join(cyc_sf_unique, cyc_data) %>%
    dplyr::select(REEF_ID, year, Hs4MW_value, geometry) %>%
    suppressMessages()
  cyc_data_sf$year <- as.numeric(cyc_data_sf$year)
  
  dhw_data_sf <- left_join(dhw_sf_unique, dhw_data) %>%
    dplyr::select(REEF_ID, year, annMaxDHW_value, geometry) %>%
    suppressMessages()
  dhw_data_sf$year <- as.numeric(dhw_data_sf$year)
  
  # Intersect cots data and aims data
  ## Returns list of distinct aims geometry points and index in unique 
  ## cots disturbance polygon dataset if within the bounds
  distinct_sample_geom <- distinct(df, geometry, .keep_all = TRUE)
  all_reefs_sf_cots <- st_intersects(distinct_sample_geom,  
                                     cots_sf_unique)
  
  all_reefs_sf_cyc <- st_intersects(distinct_sample_geom,  
                                    cyc_sf_unique)
  
  all_reefs_sf_dhw <- st_intersects(distinct_sample_geom,  
                                    dhw_sf_unique)
  
  # Find where there are no geom matches for distinct aims geometry points
  ## Don't need to do this for cyc and dhw as they result in the same indices
  sample_sf_NA_index <- which(lengths(all_reefs_sf_cots) == 0)
  
  # Find nearest cots obs for NA aims geometry points
  ## Returns the index of the nearest neighbour in cots_sf_unique$geometry
  nearest_val <- st_nn(distinct_sample_geom$geometry[sample_sf_NA_index], 
                       cots_sf_unique$geometry,
                       k = 1,
                       progress = FALSE) %>% suppressMessages()
  
  # Put these nearest neighbours in all_reefs_sf_cots
  all_reefs_sf_cots[sample_sf_NA_index] <- nearest_val
  all_reefs_sf_cyc[sample_sf_NA_index] <- nearest_val
  all_reefs_sf_dhw[sample_sf_NA_index] <- nearest_val
  
  # Join all_reefs_sf_cots to distinct_aims_geom
  geom_closest_cots <- distinct_sample_geom %>% 
    mutate(closest_index = all_reefs_sf_cots)
  geom_closest_cyc <- distinct_sample_geom %>% 
    mutate(closest_index = all_reefs_sf_cyc)
  geom_closest_dhw <- distinct_sample_geom %>% 
    mutate(closest_index = all_reefs_sf_dhw)
  
  # Loop through each row in all_reefs_sf
  for (i in seq_len(nrow(df))) {
    ## COTS ##
    # Get the index of the nearest neighbor in cots_sf_unique
    nearest_index <- geom_closest_cots$closest_index[geom_closest_cots$point_id == df$point_id[i]] %>% 
      unlist()
    
    # Get the REEF_ID of the nearest neighbor
    nearest_point_id <- cots_sf_unique$REEF_ID[nearest_index]
    
    # Find the corresponding rows in cots_data_sf based on REEF_ID and year
    cots_row <- cots_data_sf %>%
      filter(REEF_ID == nearest_point_id)
    
    ## WAVE HEIGHT / CYCLONES ##
    # Get the index of the nearest neighbor in cots_sf_unique
    nearest_index <- geom_closest_cyc$closest_index[geom_closest_cyc$point_id == df$point_id[i]] %>% 
      unlist()
    
    # Get the REEF_ID of the nearest neighbor
    nearest_point_id <- cyc_sf_unique$REEF_ID[nearest_index]
    
    # Find the corresponding rows in cyc_data_sf based on REEF_ID and year
    cyc_row <- cyc_data_sf %>%
      filter(REEF_ID == nearest_point_id)
    
    ## DEGREE HEATING WEEKS / HEATWAVES ##
    # Get the index of the nearest neighbor in cots_sf_unique
    nearest_index <- geom_closest_dhw$closest_index[geom_closest_dhw$point_id == df$point_id[i]]%>% 
      unlist()
    
    # Get the REEF_ID of the nearest neighbor
    nearest_point_id <- dhw_sf_unique$REEF_ID[nearest_index]
    
    # Find the corresponding rows in dhw_data_sf based on REEF_ID and year
    dhw_row <- dhw_data_sf %>%
      filter(REEF_ID == nearest_point_id)
    
    dists <- cots_row %>%
      mutate(Hs4MW_value = cyc_row$Hs4MW_value,
             annMaxDHW_value = dhw_row$annMaxDHW_value)
    dists <- dists[dists$year < 2017 ,c("REEF_ID", "geometry", "year", 
                                        "COTS_value", "Hs4MW_value", 
                                        "annMaxDHW_value")]
    
    # Add column for whether there's a disturbance that year
    dists <- dists %>%
      mutate(isDisturbed = (dists$COTS_value >= cots_dist | 
                              dists$Hs4MW_value >= cyc_dist |
                              dists$annMaxDHW_value >= dhw_dist))
    
    if (isTimeBased) {
      recovery_time <- recov_yrs
    } else {
      # Get recovery rate at reef from sample
      # recovery time is 1/rate
      recovery_time <- 1/as.numeric(df$r_single[i])
    }
    
    skipNext <- 0
    single_dists <- 0
    compound_dists <- 0
    for (row in which(dists$isDisturbed)) {
      if (skipNext > 0) {
        skipNext <- skipNext - 1
      } else {
        # Initialize variables to track compound disturbances
        compound_count <- 0
        current_year <- dists$year[row]
        
        # Check for disturbances within the recovery_time
        while (any(dists$isDisturbed[dists$year > current_year & 
                                     dists$year <= current_year + recovery_time])) {
          temp_dists <- dists[dists$year > current_year & 
                                dists$year <= current_year + recovery_time,]
          
          # Increment compound_count and update current_year
          compound_count <- compound_count + nrow(temp_dists)
          current_year <- max(temp_dists$year)
        }
        
        if (compound_count > 0) {
          compound_dists <- compound_dists + 1
          skipNext <- compound_count - 1  # Subtract 1 to account for the current disturbance
        } else {
          single_dists <- single_dists + 1
        }
      }
    }
    
    years_obsvd <- dists$year[nrow(dists)] - dists$year[1]
    
    df$d_single[i] <- single_dists/years_obsvd
    df$d_comp[i] <- compound_dists/years_obsvd
    
  } 
  
  # Return 
  return(df)
}

