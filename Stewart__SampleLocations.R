# Ideally this script will create 1000 samples, find the nearest environmental obs
# translate them to disturbance obs, find the nearest coral obs, calculate 
# d_{i,k} and r_{i,k} and save as a data frame to load into the main program

########### SET UP THE WORKSPACE ###########
# Clear environment
rm(list = ls())

# Clear plots
if(!is.null(dev.list())) dev.off()

# Clear commands
cat("\014")

############################################

##### LOAD LIBRARIES FUNCTIONS & DATA ######
start_time <- Sys.time()
# Load Libraries
library(sf)               # spatial feature handling
library(lwgeom)           # spatial feature handling
library(dplyr)            # dataframe manipulation
library(varhandle)        # handles variables
library(ompr)             # Optimization Modeling Package for R
library(ompr.roi)         # Optimization Modeling Package for R, ROI solver
library(reshape2)         # reshape data
library(nngeo)            # k-Nearest Neighbor Join for Spatial Data
library(lubridate)        # easy and fast parsing of date-times
library(tidyr)            # tidy data
library(ggpubr)           # arranging plots
library(ROI.plugin.glpk)  # GNU Linear Programming Kit
library(ROI)              # R Optimization Infrastructure
library(ggplot2)          # creates plots
library(latex2exp)        # LaTeX for ggplot
library(parallel)

# Set the mphil path 
mphil_path <- "../OneDrive - Queensland University of Technology/Documents/MPhil"

# Set the data path 
data_path <- paste0(mphil_path, "/Data")

# Load disturbances and coral obs (.rds made in code/obs_dist_combine.R)
all_reefs_sf <- readRDS(file = paste0(data_path,
                                      "/all_reefs_sf_gaps_filled.rds"))

# Load the shapefile
shapefile_path <- paste0(data_path,
                         "/GBRMPA/Management_Areas_of_the_Great_Barrier",
                         "_Reef_Marine_Park/Management_Areas_of_the_Great_",
                         "Barrier_Reef_Marine_Park.shp")
sector_boundaries <- st_read(shapefile_path,
                             quiet = TRUE)

# Import Environment Data
cots_data_raw <- read.csv(paste0(data_path,
                                 "/GBRdata/CoTS_data.csv"))
cyc_data_raw <- read.csv(paste0(data_path,
                                "/GBRdata/Cyclones_data.csv"))
dhw_data_raw <- read.csv(paste0(data_path,
                                "/GBRdata/DHW_data.csv"))
############################################

###### SET VARIABLES TO CHANGE ON RUN ######
num_samples <- 100000
n_threads <- detectCores()
reefs_per_mgmt <- num_samples/nrow(sector_boundaries)
############################################

####### INITIALISE VARIABLES AND DFS #######
# Minimum values considered a "disturbance" in a year
cots_dist <- 1 # cots per tow
cyc_dist <- 10 # hours of 4m wave height
dhw_dist <- 4  # degree heating weeks
# Initialise a list of random points
random_points <- data.frame(point_id = seq_len(num_samples),
                            Location = integer(length = num_samples),
                            Sector = rep(sector_boundaries$AREA_DESCR, 
                                         each = reefs_per_mgmt))
############################################

########## SAMPLE REEF LOCATIONS ###########
# Get random locations in each management area for sample reefs 
for (row in 1:nrow(sector_boundaries)) {
  from_row <- (row - 1)*reefs_per_mgmt + 1
  to_row <- row*reefs_per_mgmt
  random_points$Location[from_row:to_row] <- sector_boundaries$geometry[row] %>%
    st_sample(size = reefs_per_mgmt, 
              type = "random")
}
############################################

########### GET NEAREST ENV OBS ############
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
distinct_sample_geom <- distinct(random_points, Location, .keep_all = TRUE) %>%
  st_as_sf(crs = st_crs(sector_boundaries))

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
nearest_val <- st_nn(distinct_sample_geom$Location[sample_sf_NA_index], 
                     cots_sf_unique$geometry,
                     k = 1,
                     progress = FALSE,
                     parallel = n_threads) %>% suppressMessages()

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

# Initialise an empty data frame
columns <- c("point_id", "point_loc", "REEF_ID",
             "year", "COTS_value", "Hs4MW_value", 
             "annMaxDHW_value", "isDisturbed")
random_reefs_df <- matrix(nrow = 0, 
                          ncol = length(columns)) %>%
  data.frame()
colnames(random_reefs_df) <- columns

# Loop through each row in all_reefs_sf
for (i in seq_len(nrow(random_points))) {
  ## COTS ##
  # Get the index of the nearest neighbor in cots_sf_unique
  nearest_index <- geom_closest_cots$closest_index[geom_closest_cots$point_id == random_points$point_id[i]] %>% 
    unlist()
  
  # Get the REEF_ID of the nearest neighbor
  nearest_point_id <- cots_sf_unique$REEF_ID[nearest_index]
  
  # Find the corresponding rows in cots_data_sf based on REEF_ID and year
  cots_row <- cots_data_sf %>%
    filter(REEF_ID == nearest_point_id)
  
  ## WAVE HEIGHT / CYCLONES ##
  # Get the index of the nearest neighbor in cots_sf_unique
  nearest_index <- geom_closest_cyc$closest_index[geom_closest_cyc$point_id == random_points$point_id[i]] %>% 
    unlist()
  
  # Get the REEF_ID of the nearest neighbor
  nearest_point_id <- cyc_sf_unique$REEF_ID[nearest_index]
  
  # Find the corresponding rows in cyc_data_sf based on REEF_ID and year
  cyc_row <- cyc_data_sf %>%
    filter(REEF_ID == nearest_point_id)
  
  ## DEGREE HEATING WEEKS / HEATWAVES ##
  # Get the index of the nearest neighbor in cots_sf_unique
  nearest_index <- geom_closest_dhw$closest_index[geom_closest_dhw$point_id == random_points$point_id[i]]%>% 
    unlist()
  
  # Get the REEF_ID of the nearest neighbor
  nearest_point_id <- dhw_sf_unique$REEF_ID[nearest_index]
  
  # Find the corresponding rows in dhw_data_sf based on REEF_ID and year
  dhw_row <- dhw_data_sf %>%
    filter(REEF_ID == nearest_point_id)
  dists <- cots_row %>%
    mutate(Hs4MW_value = cyc_row$Hs4MW_value,
           annMaxDHW_value = dhw_row$annMaxDHW_value,
           point_id = rep(i, nrow(cots_row)),
           point_loc = rep(random_points$Location[i], 
                           nrow(cots_row)))
  dists <- dists %>%
    mutate(isDisturbed = (dists$COTS_value >= cots_dist | 
                            dists$Hs4MW_value >= cyc_dist |
                            dists$annMaxDHW_value >= dhw_dist)) %>%
    filter(year < 2017) %>%
    as.data.frame() %>%
    subset(select = columns)
  
  # Add to big data frame
  from_row_num <- nrow(random_reefs_df) + 1
  to_row_num <- nrow(random_reefs_df) + nrow(dists)
  random_reefs_df[from_row_num:to_row_num,] <- dists
}

############################################

################ SAVE DATA #################
# Save to an rds file
saveRDS(random_reefs_df, file = paste0(data_path,
                                       "/random_reefs_df.rds"))
end_time <- Sys.time()
end_time - start_time
############################################
