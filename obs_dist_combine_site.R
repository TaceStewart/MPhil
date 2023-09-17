## Clear environment
rm(list = ls())

## Clear commands
cat("\014")

# Set the working directory to source location
setwd("~/MPhil/Code")

library(nngeo)
library(lubridate)
library(tidyr)
library(dplyr)           

# Load the shapefile
shapefile_path <- "../Data/GBRMPA/Management_Areas_of_the_Great_Barrier_Reef_Marine_Park/Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.shp"
sector_boundaries <- st_read(shapefile_path,
                             quiet = TRUE)

## Import Data ##
cots_data_raw <- read.csv("../Data/GBRdata/CoTS_data.csv")
cyc_data_raw <- read.csv("../Data/GBRdata/Cyclones_data.csv")
dhw_data_raw <- read.csv("../Data/GBRdata/DHW_data.csv")

## Import AIMS by-site data ##
aims_site <- read.csv("../Data/AIMS_manta-tow-by-reef/ltmp_data_by_site_OA.csv")

## Import MMP Data by site ##
mmp_data <- read.csv("../Data/AIMS_MMP/mmp_hc_sc_a_by_site.csv") 

## Format reef IDs and names, add report year column
mmp_data$MMP_SITE_NAME <- gsub("Rf", "Reef", mmp_data$MMP_SITE_NAME)
mmp_data <- mutate(mmp_data,
                   REEF_SITE_NAME = paste(MMP_SITE_NAME, "Site", SITE_NO),
                   REPORT_YEAR = year(SAMPLE_DATE),
                   MONTH = month(SAMPLE_DATE))

## Multiply aims site cover by 100 as it is in decimal form
aims_site$COVER <- aims_site$COVER * 100

## Calculate report year in mmp to be comparable to aims
for (row in 1:nrow(mmp_data)) {
  if (mmp_data$MONTH[row] > 6) {
    mmp_data$REPORT_YEAR[row] <- mmp_data$REPORT_YEAR[row] + 1
  }
}

## Combine LTMP and MMP (select only Hard Cover from MMP)
hc_indices <- mmp_data$GROUP_CODE == "Hard Coral"
d5_indices <- mmp_data$DEPTH == 5
all_reefs <- data.frame(REEF_NAME = c(aims_site$REEF, mmp_data$REEF_SITE_NAME[hc_indices & d5_indices]),
                        DEPTH = c(rep(NA, nrow(aims_site)), mmp_data$DEPTH[hc_indices & d5_indices]),
                        LATITUDE = c(aims_site$LATITUDE, mmp_data$LATITUDE[hc_indices & d5_indices]),
                        LONGITUDE = c(aims_site$LONGITUDE, mmp_data$LONGITUDE[hc_indices & d5_indices]),
                        YEAR = c(aims_site$REPORT_YEAR, year(mmp_data$SAMPLE_DATE[hc_indices & d5_indices])),
                        COVER = c(aims_site$COVER, mmp_data$COVER[hc_indices & d5_indices]),
                        PROGRAM = c(rep("LTMP", nrow(aims_site)), rep("MMP", sum(hc_indices & d5_indices))))

# Create spatial points for the reefs using the longitude and latitude coordinates
all_reefs_sf <- st_as_sf(all_reefs, 
                         coords = c("LONGITUDE", "LATITUDE"), 
                         crs = st_crs(sector_boundaries))

# Perform a spatial join to determine which area each reef belongs to
all_reefs_sf <- st_join(all_reefs_sf, sector_boundaries[, "AREA_DESCR"])

# Remove any coral obs after 2016 as not all disturbance data covers it
all_reefs_sf <- all_reefs_sf[all_reefs_sf$YEAR < 2017,]

# Transform columns to rows, summarise cots at each Reef ID
cots_data <- cots_data_raw %>%
  pivot_longer(cols = starts_with("COTS"),
               names_to = "year",
               values_to = "COTS_value",
               values_drop_na = FALSE) %>%
  mutate(year = stringr::str_extract(year, "\\d{4}")) %>%
  group_by(REEF_ID,year) %>%
  summarize(COTS_value = sum(COTS_value, na.rm = TRUE))

cyc_data <- cyc_data_raw %>%
  pivot_longer(cols = starts_with("Hs4MW"),
               names_to = "year",
               values_to = "Hs4MW_value",
               values_drop_na = FALSE) %>%
  mutate(year = stringr::str_extract(year, "\\d{4}")) %>%
  group_by(REEF_ID,year) %>%
  summarize(Hs4MW_value = mean(Hs4MW_value, na.rm = TRUE))

dhw_data <- dhw_data_raw %>%
  pivot_longer(cols = starts_with("annMaxDHW"),
               names_to = "year",
               values_to = "annMaxDHW_value",
               values_drop_na = FALSE) %>%
  mutate(year = stringr::str_extract(year, "\\d{4}")) %>%
  group_by(REEF_ID,year) %>%
  summarize(annMaxDHW_value = mean(annMaxDHW_value, na.rm = TRUE))

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
  dplyr::select(REEF_ID, year, COTS_value, geometry)
cots_data_sf$year <- as.numeric(cots_data_sf$year)

cyc_data_sf <- left_join(cyc_sf_unique, cyc_data) %>%
  dplyr::select(REEF_ID, year, Hs4MW_value, geometry)
cyc_data_sf$year <- as.numeric(cyc_data_sf$year)

dhw_data_sf <- left_join(dhw_sf_unique, dhw_data) %>%
  dplyr::select(REEF_ID, year, annMaxDHW_value, geometry)
dhw_data_sf$year <- as.numeric(dhw_data_sf$year)

# Intersect cots data and aims data
## Returns list of distinct aims geometry points and index in unique 
## cots disturbance polygon dataset if within the bounds
distinct_aims_geom <- distinct(all_reefs_sf, geometry, .keep_all = TRUE)
all_reefs_sf_cots <- st_intersects(distinct_aims_geom,  
                                   cots_sf_unique)

all_reefs_sf_cyc <- st_intersects(distinct_aims_geom,  
                                  cyc_sf_unique)

all_reefs_sf_dhw <- st_intersects(distinct_aims_geom,  
                                  dhw_sf_unique)

# Find where there are no geom matches for distinct aims geometry points
## Don't need to do this for cyc and dhw as they result in the same indices
AIMS_sf_NA_index <- which(lengths(all_reefs_sf_cots) == 0)

# Find nearest cots obs for NA aims geometry points
## Returns the index of the nearest neighbour in cots_sf_unique$geometry
nearest_val <- st_nn(distinct_aims_geom$geometry[AIMS_sf_NA_index], 
                     cots_sf_unique$geometry,
                     k = 1)

# Put these nearest neighbours in all_reefs_sf_cots
all_reefs_sf_cots[AIMS_sf_NA_index] <- nearest_val
all_reefs_sf_cyc[AIMS_sf_NA_index] <- nearest_val
all_reefs_sf_dhw[AIMS_sf_NA_index] <- nearest_val

# Join all_reefs_sf_cots to distinct_aims_geom
geom_closest_cots <- distinct_aims_geom %>% 
  mutate(closest_index = all_reefs_sf_cots)
geom_closest_cyc <- distinct_aims_geom %>% 
  mutate(closest_index = all_reefs_sf_cyc)
geom_closest_dhw <- distinct_aims_geom %>% 
  mutate(closest_index = all_reefs_sf_dhw)

# Initialize an empty vector to store the COTS values
all_reefs_sf$COTS_value <- NA
all_reefs_sf$Hs4MW_value <- NA
all_reefs_sf$DHW_value <- NA

# Loop through each row in all_reefs_sf
for (i in seq_len(nrow(all_reefs_sf))) {
  ## COTS ##
  # Get the index of the nearest neighbor in cots_sf_unique
  nearest_index <- unlist(geom_closest_cots$closest_index[geom_closest_cots$REEF_NAME == all_reefs_sf$REEF_NAME[i]])
  
  # Get the REEF_ID of the nearest neighbor
  nearest_reef_id <- cots_sf_unique$REEF_ID[nearest_index]
  
  # Find the corresponding rows in cots_data_sf based on REEF_ID and year
  cots_row <- cots_data_sf %>%
    filter(REEF_ID == nearest_reef_id)
  
  ## WAVE HEIGHT / CYCLONES ##
  # Get the index of the nearest neighbor in cyc_sf_unique
  nearest_index <- unlist(geom_closest_cyc$closest_index[geom_closest_cyc$REEF_NAME == all_reefs_sf$REEF_NAME[i]])
  
  # Get the REEF_ID of the nearest neighbor
  nearest_reef_id <- cyc_sf_unique$REEF_ID[nearest_index]
  
  # Find the corresponding rows in cyc_data_sf based on REEF_ID and year
  cyc_row <- cyc_data_sf %>%
    filter(REEF_ID == nearest_reef_id)
  
  ## DEGREE HEATING WEEKS / HEATWAVES ##
  # Get the index of the nearest neighbor in cots_sf_unique
  nearest_index <- unlist(geom_closest_dhw$closest_index[geom_closest_dhw$REEF_NAME == all_reefs_sf$REEF_NAME[i]])
  
  # Get the REEF_ID of the nearest neighbor
  nearest_reef_id <- dhw_sf_unique$REEF_ID[nearest_index]
  
  # Find the corresponding rows in dhw_data_sf based on REEF_ID and year
  dhw_row <- dhw_data_sf %>%
    filter(REEF_ID == nearest_reef_id)
  # TO DO: MAKE THIS WORK FOR SAME SITE PREV YEAR
  # If this is not the first observation in a reef
  if (i != 1 && all_reefs_sf$REEF_NAME[i-1] == all_reefs_sf$REEF_NAME[i]) {
    # Get the years for the last obs until the current row in all_reefs_sf
    year_current_row <- all_reefs_sf$YEAR[i-1:i]
    
    # Find the matching COTS observation for the current reef and years
    matching_cots_value <- max(cots_row$COTS_value[cots_row$year %in% 
                                                     year_current_row],
                               na.rm = TRUE)
    
    # Find the matching cyc observation for the current reef and years
    matching_cyc_value <- max(cyc_row$Hs4MW_value[cyc_row$year %in% 
                                                    year_current_row],
                              na.rm = TRUE)
    
    # Find the matching dhw observation for the current reef and year
    matching_dhw_value <- max(dhw_row$annMaxDHW_value[dhw_row$year %in% 
                                                        year_current_row],
                              na.rm = TRUE)
  } else {
    # Get the current year of obs
    year_current_row <- all_reefs_sf$YEAR[i]
    
    # Find the matching COTS observation for the current reef and year
    matching_cots_value <- cots_row$COTS_value[cots_row$year == year_current_row]
    
    # Find the matching cyc observation for the current reef and years
    matching_cyc_value <- cyc_row$Hs4MW_value[cyc_row$year == year_current_row]
    
    # Find the matching dhw observation for the current reef and year
    matching_dhw_value <- dhw_row$annMaxDHW_value[dhw_row$year == year_current_row]
  }
  
  if (is.na(matching_cots_value)) {
    temp <- cots_data_sf[cots_data_sf$year == year_current_row,]
    temp <- na.omit(temp)
    
    # Find the nearest not-NA value
    nearest_val <- unlist(st_nn(all_reefs_sf$geometry[i], 
                                temp$geometry,
                                k = 1,
                                progress = FALSE))
    
    # Find the matching COTS observation for the current reef and year
    matching_cots_value <- temp$COTS_value[nearest_val]
  }
  
  if (is.na(matching_cyc_value)) {
    temp <- cyc_data_sf[cyc_data_sf$year == year_current_row,]
    temp <- na.omit(temp)
    
    if (nrow(temp) > 0) {
      # Find the nearest not-NA value
      nearest_val <- unlist(st_nn(all_reefs_sf$geometry[i], 
                                  temp$geometry,
                                  k = 1,
                                  progress = FALSE))
      
      # Find the matching cyc observation for the current reef and year
      matching_cyc_value <- temp$Hs4MW_value[nearest_val]
    } else {
      matching_cyc_value <- NA
    }
  }
  
  if (is.na(matching_dhw_value)) {
    temp <- dhw_data_sf[dhw_data_sf$year == year_current_row,]
    temp <- na.omit(temp)
    
    if (nrow(temp) > 0) {
      # Find the nearest not-NA value
      nearest_val <- unlist(st_nn(all_reefs_sf$geometry[i], 
                                  temp$geometry,
                                  k = 1,
                                  progress = FALSE))
      
      # Find the matching dhw observation for the current reef and year
      matching_dhw_value <- temp$DHW_value[nearest_val]
    } else {
      matching_dhw_value <- NA
    }
  }
  
  # Update the COTS_value column for the matched year and reef
  all_reefs_sf$COTS_value[i] <- matching_cots_value
  
  # Update the cyc_value column for the matched year and reef
  all_reefs_sf$Hs4MW_value[i] <- matching_cyc_value
  
  # Update the COTS_value column for the matched year and reef
  all_reefs_sf$DHW_value[i] <- matching_dhw_value
} 

# Save to an rds file
saveRDS(all_reefs_sf, file = "../Data/all_reef_sites_sf.rds")