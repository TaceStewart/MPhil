# Clear environment
rm(list = ls())

# Clear commands
cat("\014")

# Load libraries
library(nngeo)
library(lubridate)
library(tidyr)
library(dplyr) 

# Set the data path 
data_path <- "../OneDrive - Queensland University of Technology/Documents/MPhil/Data"

# figure out where the gaps are in ltmp dataset summd by reef, that can be 
# filled by the by-site
all_reefs <- readRDS(file = paste0(data_path,
                                   "/all_reefs_sf.rds"))
all_sites <- readRDS(file = paste0(data_path,
                                   "/all_reef_sites_sf.rds"))

ltmp_reefs <- all_reefs[all_reefs$PROGRAM == "LTMP",]
ltmp_sites <- all_sites[all_sites$PROGRAM == "LTMP",]

ltmp_sites$COVER <- ltmp_sites$COVER

unique_reefs <- unique(ltmp_reefs$REEF_NAME)
unique_sites <- ltmp_sites$REEF_NAME %>% unique() 

unique_reef_locs <- distinct(ltmp_reefs, geometry)
unique_site_locs <- distinct(ltmp_sites, geometry)

nearest_val <- st_nn(ltmp_reefs$geometry, 
                     ltmp_sites$geometry,
                     k = 1,
                     maxdist = 2500)
nearest_indices <- which(nearest_val > 0)

ltmp_reefs$nearest_val <- NA
ltmp_reefs$nearest_val[nearest_indices] <- unlist(nearest_val)

ltmp_sites$REEF_NAME <- gsub(" Site \\d{1}", "", ltmp_sites$REEF_NAME)

count <- 0
reef_names <- unique(ltmp_reefs$REEF_NAME[!is.na(ltmp_reefs$nearest_val)])
for (reef_name in reef_names) {
  reef_obs <- ltmp_reefs[ltmp_reefs$REEF_NAME == reef_name,]
  site_name <- ltmp_sites$REEF_NAME[unlist(reef_obs$nearest_val[1])]
  site_obs <- ltmp_sites[ltmp_sites$REEF_NAME == site_name,]
  
  reef_years <- unique(reef_obs$YEAR)
  site_years <- unique(site_obs$YEAR)
  
  if (any(!(site_years %in% reef_years))) {
    count <- count + sum(!(site_years %in% reef_years))
    
    site_means <- site_obs %>%
      group_by(YEAR) %>%
      summarize(meanCOVER = mean(COVER),
                meanCOTS = mean(COTS_value),
                meanWAVE = mean(Hs4MW_value),
                meanDHW = mean(DHW_value)) %>%
      ungroup()
    
    # reefplot <- ggplot() +
    #   geom_point(data = reef_obs,
    #              mapping = aes(x = YEAR, y = COVER),
    #              colour = "red") +
    #   geom_line(data = reef_obs,
    #             mapping = aes(x = YEAR, y = COVER),
    #             colour = "red") +
    #   geom_point(data = site_means,
    #              mapping = aes(x = YEAR, y = meanCOVER),
    #              colour = "blue") +
    #   geom_line(data = site_means,
    #             mapping = aes(x = YEAR, y = meanCOVER),
    #             colour = "blue") +
    #   scale_y_continuous(limits = c(0,100))
    
    for (index in which(!(site_years %in% reef_years))) {
      all_reefs[nrow(all_reefs) + 1,] <- data.frame(reef_name, 
                                                    site_obs$DEPTH[1],
                                                    site_means$YEAR[index],
                                                    site_means$meanCOVER[index],
                                                    site_obs$PROGRAM[index],
                                                    site_obs$AREA_DESCR[index],
                                                    reef_obs$geometry[1],
                                                    site_means$meanCOTS[index],
                                                    site_means$meanWAVE[index],
                                                    site_means$meanDHW[index])
    }
  }
}

# Sort it by reef then by year
all_reefs_sorted <- all_reefs[order(all_reefs$REEF_NAME, all_reefs$YEAR),]

# Save to an rds file
saveRDS(all_reefs_sorted, file = paste0(data_path,
                                        "/all_reefs_sf_gaps_filled.rds"))