########### LOAD LIBRARIES ######
library(leaflet)
library(maps)
library(raster)
library(tidyverse)
library(ggspatial)
library(dataaimsr)
library(gisaimsr)
library(ggrepel)
library(sf)       # for spatial data handling
library(ggplot2)  # for data visualization
library(leaflet.providers)
library(htmlwidgets)
library(htmltools)
library(webshot)
library(hexbin)
##########################################

############## LOAD DATA #########

# Set the mphil path 
mphil_path <- "../OneDrive - Queensland University of Technology/Documents/MPhil"

# Set the data path 
data_path <- paste0(mphil_path, "/Data")

# Set the figure output path 
out_path <- paste0(mphil_path, "/Figures/Sensitivity Outputs")

# Load the shapefile
shapefile_path <- paste0(data_path, 
                         "/GBRMPA/Management_Areas_of_the_Great_Barrier_Reef_Marine_Park/Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.shp")
map_data <- st_read(shapefile_path)

# Load .RData
#load("~/MPhil/Code/Parameter Run Environments/RecovBased0.75th0.1mgmt.RData")
load(paste0(data_path, 
            "/Parameter Run Environments/TimeBased5yrs0.1mgmt.RData"))
#load("Parameter Run Environments/RecovBased0.75th0.1mgmt.RData")
##########################################

#### FOR THESIS DATA CHAPTER ####
allreefs_xy <- st_coordinates(all_reefs_sf)
all_reefs_sf$X <- allreefs_xy[,1]
all_reefs_sf$Y <- allreefs_xy[,2]
landscape_dims <- c(8,4)
portrait_dims <- c(5,6)

##### Management Areas (Crop after) ######
rebe_coords_xy <- all_reefs_sf$geometry[which(grepl("Rebe", 
                                                    all_reefs_sf$REEF_NAME,
                                                    ignore.case = T))] %>% 
  st_coordinates()
coordinates_xy <- st_coordinates(all_reefs_sf$geometry)
map_data$NUM_REEFS <- 0
unique_rfs <- unique(all_reefs_sf[, c("REEF_NAME", "AREA_DESCR")])
for (sector_row in 1:nrow(map_data)) {
  map_data$NUM_REEFS[sector_row] <- sum(unique_rfs$AREA_DESCR == 
                                          map_data$AREA_DESCR[sector_row])
}
leaflet() %>%
  setView(lng = mean(coordinates_xy[,1]), 
          lat = min(coordinates_xy[,2]) +5, 
          zoom = 5) %>%
  addProviderTiles("Esri.WorldGrayCanvas") %>%
  addPolygons(data = map_data, 
              label = ~paste0(gsub(" Management Area", 
                                   "",
                                   AREA_DESCR), " (n = ", NUM_REEFS, ")"),
              labelOptions = labelOptions(noHide = T, 
                                          textOnly = TRUE, 
                                          direction = "right",
                                          style = list(
                                            "color" = "grey",
                                            "font-weight" = "bold",
                                            "font-size" = "14px"),
                                          offset = c(38, -15)),
              color = "grey",
              opacity = 1,
              weight = 2) %>%
  addLabelOnlyMarkers(lng = 140,
                      lat = -25,
                      label = "Queensland",
                      labelOptions = labelOptions(noHide = T, 
                                                  textOnly = TRUE,
                                                  style = list(
                                                    "color" = "grey",
                                                    "font-weight" = "bold",
                                                    "font-size" = "12px"))) %>%
  addCircleMarkers(lng = mean(rebe_coords_xy[,1]),
                   lat = mean(rebe_coords_xy[,2]),
                   color = "lightcoral",
                   opacity = 1,
                   weight = 4,
                   radius = 5,
                   fill = FALSE) %>%
  addLabelOnlyMarkers(lng = mean(rebe_coords_xy[,1] + 1),
                      lat = mean(rebe_coords_xy[,2]),
                      label = "Rebe Reef",
                      labelOptions = labelOptions(noHide = T, 
                                                  textOnly = TRUE, 
                                                  direction = "right",
                                                  style = list(
                                                    "color" = "lightcoral",
                                                    "font-weight" = "bold",
                                                    "font-size" = "18px"))) %>%
  saveWidget("gbrmpaMA.html")
webshot("gbrmpaMA.html", 
        "../MPhil Thesis/Figures/Data Chapter/Aus_GBR_FourMAs.png",
        vwidth = 700, vheight = 500, zoom = 8)

##########################################

##### Aims monitoring locs by program ####
unique_rfs <- unique(all_reefs_sf[, c("REEF_NAME", "PROGRAM", "geometry")])
ggplot() +
  geom_sf(data = map_data, 
          lwd = 0.75,
          color = "darkgrey") +
  theme_classic() +
  geom_sf(data = unique_rfs,
          pch = 19,
          mapping = aes(colour = PROGRAM),
          alpha = 0.75) +
  labs(x = "Longitude",
       y = "Latitude",
       colour = "AIMS Program") + 
  theme(legend.position = c(0.05, 0.05),
        legend.box.background = element_rect(linewidth=1, colour = "darkgrey"),
        legend.justification = c("left", "bottom"))
ggsave("../MPhil Thesis/Figures/Data Chapter/AIMSLocByProgram.png", 
       plot = last_plot(), width=portrait_dims[1], height=portrait_dims[2])
##########################################

###### Box plot of cover by program ######
mean_trend <- all_reefs_sf %>%
  group_by(YEAR) %>%
  summarise(mean = mean(COVER))
ggplot() +
  theme_classic() +
  geom_boxplot(data = all_reefs_sf,
               mapping = aes(PROGRAM, 
                             COVER),
               linewidth = 1,
               color = "azure4") +
  labs(x = "Program",
       y = "Cover (%)") +
  scale_y_continuous(limits = c(0,100))
ggsave("../MPhil Thesis/Figures/Data Chapter/CoralCoverBoxPlotByProgram.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

########### Single dist spatial ##########
side_size <- 0.5
ggplot() +
  theme_classic() +
  geom_sf(data = map_data,
          size = 0.25,
          color = "azure4",
          fill = NA) +
  geom_hex(data = all_reefs_sf[all_reefs_sf$single_or_compound == "Single",],
           aes(X, Y),
           binwidth = c(side_size, side_size),
           na.rm = TRUE) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Number of \nSingle \nDisturbances")
if (isTimeBased) {
  ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/GBRSingleDistsTimeBased", 
                recov_yrs, "yr", mgmt_benefit, "mgmt.png"),
         plot = last_plot(), width=5, height=5)
} else {
  ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/GBRSingleDistsRecovBased", 
                recov_th, "th", mgmt_benefit, "mgmt.png"),
         plot = last_plot(), width=5, height=5)
}
##########################################

########## Compound dist spatial #########
side_size <- 0.5
ggplot() +
  theme_classic() +
  geom_sf(data = map_data,
          size = 0.25,
          color = "azure4",
          fill = NA) +
  geom_hex(data = all_reefs_sf[all_reefs_sf$single_or_compound == "Compound",],
           aes(X, Y),
           binwidth = c(side_size, side_size),
           na.rm = TRUE) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Number of \nCompound \nDisturbances")
if (isTimeBased) {
  ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/GBRCompoundDistsTimeBased", 
                recov_yrs, "yr", mgmt_benefit, "mgmt.png"),
         plot = last_plot(), width=5, height=5)
} else {
  ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/GBRCompoundDistsRecovBased", 
                recov_th, "th", mgmt_benefit, "mgmt.png"),
         plot = last_plot(), width=5, height=5)
}
summarised_reefs <- all_reefs_sf %>%
  group_by(REEF_NAME) %>%
  summarise(total = n())
##########################################

##### Rebe Reef location (Crop after) #####

leaflet() %>%
  setView(lng = mean(rebe_coords_xy[,1]) - 1, 
          lat = mean(rebe_coords_xy[,2]) + 1, 
          zoom = 5) %>%
  addProviderTiles("Esri.WorldGrayCanvas") %>%
  addPolygons(data = map_data, 
              popup = ~AREA_DESCR,
              color = "grey",
              opacity = 1,
              weight = 1) %>% 
  addCircleMarkers(lng = mean(rebe_coords_xy[,1]),
                   lat = mean(rebe_coords_xy[,2]),
                   color = "lightcoral",
                   opacity = 1,
                   weight = 4,
                   radius = 5,
                   fill = FALSE) %>%
  addLabelOnlyMarkers(lng = mean(rebe_coords_xy[,1] + 1),
                      lat = mean(rebe_coords_xy[,2]),
                      label = "Rebe Reef",
                      labelOptions = labelOptions(noHide = T, 
                                                  textOnly = TRUE, 
                                                  direction = "right",
                                                  style = list(
                                                    "color" = "lightcoral",
                                                    "font-weight" = "bold",
                                                    "font-size" = "18px"))) %>%
  addLabelOnlyMarkers(lng = 140,
                      lat = -25,
                      label = "Queensland",
                      labelOptions = labelOptions(noHide = T, 
                                                  textOnly = TRUE,
                                                  style = list(
                                                    "color" = "grey",
                                                    "font-weight" = "bold",
                                                    "font-size" = "12px"))) %>%
  saveWidget("rebe.html")
webshot("rebe.html", 
        "../MPhil Thesis/Figures/Data Chapter/GBR_RebeReef.png",
        vwidth = 650, vheight = 500, zoom = 8)
##########################################

######### Rebe Reef linear plot ##########
ggplot(data = all_reefs_sf[which(grepl("REBE", all_reefs_sf$REEF_NAME)),],
       mapping = aes(x = YEAR, 
                     y = COVER)) +
  theme_classic() +
  geom_point(colour = "azure4",
             size = 1.5, 
             alpha = 0.75) +
  geom_line(colour = "azure4",
            linewidth = 1,
            alpha = 0.75) +
  labs(x = "Year",
       y = "Cover (%)") + 
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(limits = c(0,50), 
                     breaks = seq(0, 50, 10))
ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefCoralCover.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

####### Rebe Reef linear with dists ######
rebe_reef <- all_reefs_sf[which(grepl("REBE", all_reefs_sf$REEF_NAME)),]
rebe_reef$yearDists <- NA
distNames <- c("CoTS Outbreak", "Wind Stress", "Heat Stress")
for (row in 1: nrow(rebe_reef)) {
  distType <- distNames[c(rebe_reef$COTS_value[row] >= cots_dist, 
                          rebe_reef$Hs4MW_value[row] >= cyc_dist, 
                          rebe_reef$DHW_value[row] >= dhw_dist)]
  if (length(distType > 1)) {
    distType <- paste(distType, collapse = ", ")
    rebe_reef$yearDists[row] <- distType
  }
}
rebe_reef$distYears <- ifelse(is.na(rebe_reef$yearDists), NA, rebe_reef$YEAR)
colors <- c("Wind Stress" = "steelblue3", 
            "Heat Stress" = "lightcoral", 
            "CoTS Outbreak" = "mediumpurple",
            "CoTS Outbreak, Wind Stress, Heat Stress" = "aquamarine3")
rebe_reef$distYears <- ifelse(is.na(rebe_reef$yearDists), NA, rebe_reef$YEAR)
ggplot(data = rebe_reef,
       mapping = aes(x = YEAR, 
                     y = COVER)) +
  theme_classic() +
  geom_point(colour = "azure4",
             size = 1.5, 
             alpha = 0.75) +
  geom_line(colour = "azure4",
            linewidth = 1,
            alpha = 0.75) +
  geom_vline(mapping = aes(xintercept = rebe_reef$distYears, 
                           colour = rebe_reef$yearDists),
             na.rm = TRUE,
             linewidth = 1)  + 
  scale_colour_manual(na.translate = F,
                      values = colors) +
  labs(x = "Year",
       y = "Cover (%)",
       colour = "Disturbance Type") + 
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(limits = c(0,50), 
                     breaks = seq(0, 50, 10))
ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefCoralCover_AllDists.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

######## Rebe Env + coral linear #########
scaleFactor <- max(rebe_reef$COVER) / max(rebe_reef$Hs4MW_value / cyc_dist * 100)
ggplot(data = rebe_reef) +
  theme_classic() +
  geom_point(mapping = aes(x = YEAR,
                           y = COVER / scaleFactor),
             size = 1.5,
             alpha = 0.75,
             color = "azure4") +
  geom_line(mapping = aes(x = YEAR,
                          y = COVER / scaleFactor),
            linewidth = 1,
            alpha = 0.75,
            color = "azure4") +
  geom_point(mapping = aes(x = YEAR, 
                           y = Hs4MW_value / cyc_dist * 100,
                           color = "Wind Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = Hs4MW_value / cyc_dist * 100,
                          color = "Wind Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = DHW_value / dhw_dist * 100,
                           color = "Heat Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = DHW_value / dhw_dist * 100,
                          color = "Heat Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = COTS_value / cots_dist * 100,
                           color = "CoTS Outbreak"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = COTS_value / cots_dist * 100,
                          color = "CoTS Outbreak"),
            alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = 100,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       y = "% of Disturbance Level",
       color = "Disturbance Type") + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = "% of Disturbance Level",
    breaks = c(0, 50, 100, 150, 200, 250, 300),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*scaleFactor, name="Coral Cover (%)")) +
  theme(axis.title.y.right=element_text(color="azure4"),
        axis.text.y.right=element_text(color="azure4"))
ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefDisturbanceWithObs.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

################ Rebe Env ################
scaleFactor <- max(rebe_reef$COVER) / max(rebe_reef$Hs4MW_value / cyc_dist * 100)
ggplot(data = rebe_reef) +
  theme_classic() +
  geom_point(mapping = aes(x = YEAR, 
                           y = Hs4MW_value / cyc_dist * 100,
                           color = "Wind Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = Hs4MW_value / cyc_dist * 100,
                          color = "Wind Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = DHW_value / dhw_dist * 100,
                           color = "Heat Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = DHW_value / dhw_dist * 100,
                          color = "Heat Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = COTS_value / cots_dist * 100,
                           color = "CoTS Outbreak"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = COTS_value / cots_dist * 100,
                          color = "CoTS Outbreak"),
            alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = 100,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       y = "% of Disturbance Level",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = "% of Disturbance Level",
    breaks = c(0, 50, 100, 150, 200, 250, 300))
ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefDisturbancesOnly.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

################ Rebe WH #################
ggplot(data = rebe_reef,
       mapping = aes(x = YEAR, 
                           y = Hs4MW_value,
                           color = "Wind Stress")) +
  theme_classic() +
  geom_point(alpha = 0.75,
             size = 1.5) +
  geom_line(alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = cyc_dist,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = TeX("Hours of $\\geq$ 4m Wave Height"))
ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefWHwithdistlvl.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

################ Rebe DHW ################
ggplot(data = rebe_reef,
       mapping = aes(x = YEAR, 
                     y = DHW_value,
                     color = "Heat Stress")) +
  theme_classic() +
  geom_point(alpha = 0.75,
             size = 1.5) +
  geom_line(alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = dhw_dist,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = TeX("Number of Degree-Heating-Weeks"))
ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefDHWwithdistlvl.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

############### Rebe CoTS ################
ggplot(data = rebe_reef,
       mapping = aes(x = YEAR, 
                     y = COTS_value,
                     color = "CoTS Outbreak")) +
  theme_classic() +
  geom_point(alpha = 0.75,
             size = 1.5) +
  geom_line(alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = cots_dist,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = TeX("Number of Crown-of-Thorns Starfish per Manta-Tow"))
ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefCOTSwithdistlvl.png", 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

######## Rebe single and compound ########
rebe_reef$distYears <- ifelse(is.na(rebe_reef$yearDists), NA, rebe_reef$YEAR)
SingCompDistYears <- ifelse(is.na(rebe_reef$single_or_compound), NA, rebe_reef$YEAR)
ggplot(data = rebe_reef,
       mapping = aes(x = YEAR, 
                     y = COVER)) +
  theme_classic() +
  geom_rect(xmin = SingCompDistYears, 
            xmax = ifelse(grepl("Unknown", rebe_reef$recovYear), 
                          Inf,
                          as.numeric(rebe_reef$recovYear)), 
            ymin = -Inf, ymax = Inf,
            mapping = aes(fill = rebe_reef$single_or_compound),
            alpha = .3,
            na.rm = TRUE) + 
  scale_fill_manual(name = "Event Type",
                      values = c("orange", "yellow"),
                      na.translate = F) +
  geom_point(colour = "azure4",
             size = 1.5,
             alpha = 0.75) +
  geom_line(colour = "azure4",
            linewidth = 1,
            alpha = 0.75) +
  geom_vline(mapping = aes(xintercept = rebe_reef$distYears), 
             colour = "gray40",
             na.rm = TRUE,
             linewidth = 1) +
  labs(x = "Year",
       y = "Cover (%)",
       colour = "Disturbance Type") + 
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(limits = c(0,50), 
                     breaks = seq(0, 50, 10))
if (isTimeBased) {
  ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefCompoundandSingleTimeBased.png", 
         plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
} else {
  ggsave("../MPhil Thesis/Figures/Data Chapter/RebeReefCompoundandSingleRecovBased.png", 
         plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
}
##########################################

reef_name <- "HASTINGS REEF"

####### reef_name linear with dists ######
obs_by_reef <- all_reefs_sf[all_reefs_sf$REEF_NAME == reef_name,]
unique(obs_by_reef$distType)
obs_by_reef$yearDists <- NA
distNames <- c("CoTS Outbreak", "Wind Stress", "Heat Stress")
for (row in 1: nrow(obs_by_reef)) {
  distType <- distNames[c(obs_by_reef$COTS_value[row] >= cots_dist, 
                          obs_by_reef$Hs4MW_value[row] >= cyc_dist, 
                          obs_by_reef$DHW_value[row] >= dhw_dist)]
  if (length(distType > 1)) {
    distType <- paste(distType, collapse = ", ")
    obs_by_reef$yearDists[row] <- distType
  }
}
obs_by_reef$distYears <- ifelse(is.na(obs_by_reef$yearDists), NA, obs_by_reef$YEAR)
colors <- c("Wind Stress" = "steelblue3", 
            "Heat Stress" = "lightcoral", 
            "CoTS Outbreak" = "mediumpurple",
            "CoTS Outbreak, CoTS Outbreak" = "darkpurple")
obs_by_reef$distYears <- ifelse(is.na(obs_by_reef$yearDists), NA, obs_by_reef$YEAR)
ggplot(data = obs_by_reef,
       mapping = aes(x = YEAR, 
                     y = COVER)) +
  theme_classic() +
  geom_point(colour = "azure4",
             size = 1.5, 
             alpha = 0.75) +
  geom_line(colour = "azure4",
            linewidth = 1,
            alpha = 0.75) +
  geom_vline(mapping = aes(xintercept = obs_by_reef$distYears, 
                           colour = obs_by_reef$yearDists),
             na.rm = TRUE,
             linewidth = 1)  + 
  scale_colour_manual(na.translate = F,
                      values = colors) +
  labs(x = "Year",
       y = "Cover (%)",
       colour = "Disturbance Type") + 
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(limits = c(0,50), 
                     breaks = seq(0, 50, 10))
ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_CoralCover_WDists.png"), 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

##### reef_name env + coral linear #####
scaleFactor <- max(obs_by_reef$COVER) / max(obs_by_reef$Hs4MW_value / cyc_dist * 100)
ggplot(data = obs_by_reef) +
  theme_classic() +
  geom_point(mapping = aes(x = YEAR,
                           y = COVER / scaleFactor),
             size = 1.5,
             alpha = 0.75,
             color = "azure4") +
  geom_line(mapping = aes(x = YEAR,
                          y = COVER / scaleFactor),
            linewidth = 1,
            alpha = 0.75,
            color = "azure4") +
  geom_point(mapping = aes(x = YEAR, 
                           y = Hs4MW_value / cyc_dist * 100,
                           color = "Wind Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = Hs4MW_value / cyc_dist * 100,
                          color = "Wind Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = DHW_value / dhw_dist * 100,
                           color = "Heat Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = DHW_value / dhw_dist * 100,
                          color = "Heat Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = COTS_value / cots_dist * 100,
                           color = "CoTS Outbreak"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = COTS_value / cots_dist * 100,
                          color = "CoTS Outbreak"),
            alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = 100,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       y = "% of Disturbance Level",
       color = "Disturbance Type") + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = "% of Disturbance Level",
    breaks = c(0, 50, 100, 150, 200, 250, 300),
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*scaleFactor, name="Coral Cover (%)")) +
  theme(axis.title.y.right=element_text(color="azure4"),
        axis.text.y.right=element_text(color="azure4"))
ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_DisturbanceWithObs.png"), 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

############# reef_name env ##############
scaleFactor <- max(obs_by_reef$COVER) / max(obs_by_reef$Hs4MW_value / cyc_dist * 100)
ggplot(data = obs_by_reef) +
  theme_classic() +
  geom_point(mapping = aes(x = YEAR, 
                           y = Hs4MW_value / cyc_dist * 100,
                           color = "Wind Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = Hs4MW_value / cyc_dist * 100,
                          color = "Wind Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = DHW_value / dhw_dist * 100,
                           color = "Heat Stress"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = DHW_value / dhw_dist * 100,
                          color = "Heat Stress"),
            alpha = 0.75,
            linewidth = 1) +
  geom_point(mapping = aes(x = YEAR, 
                           y = COTS_value / cots_dist * 100,
                           color = "CoTS Outbreak"),
             alpha = 0.75,
             size = 1.5) +
  geom_line(mapping = aes(x = YEAR, 
                          y = COTS_value / cots_dist * 100,
                          color = "CoTS Outbreak"),
            alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = 100,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       y = "% of Disturbance Level",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = "% of Disturbance Level",
    breaks = c(0, 50, 100, 150, 200, 250, 300))
ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_DisturbancesOnly.png"), 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

################ Reef WH #################
ggplot(data = obs_by_reef,
       mapping = aes(x = YEAR, 
                     y = Hs4MW_value,
                     color = "Wind Stress")) +
  theme_classic() +
  geom_point(alpha = 0.75,
             size = 1.5) +
  geom_line(alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = cyc_dist,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = TeX("Hours of $\\geq$ 4m Wave Height"))
ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_WHwithdistlvl.png"), 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

################ Reef DHW ################
ggplot(data = obs_by_reef,
       mapping = aes(x = YEAR, 
                     y = DHW_value,
                     color = "Heat Stress")) +
  theme_classic() +
  geom_point(alpha = 0.75,
             size = 1.5) +
  geom_line(alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = dhw_dist,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = TeX("Number of Degree-Heating-Weeks"))
ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_DHWwithdistlvl.png"), 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

############### Reef CoTS ################
ggplot(data = obs_by_reef,
       mapping = aes(x = YEAR, 
                     y = COTS_value,
                     color = "CoTS Outbreak")) +
  theme_classic() +
  geom_point(alpha = 0.75,
             size = 1.5) +
  geom_line(alpha = 0.75,
            linewidth = 1) +
  geom_hline(yintercept = cots_dist,
             colour = "Red",
             linetype = "dashed") +
  labs(x = "Year",
       colour = "Disturbance Type") + 
  scale_colour_manual(values = colors) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(# Features of the first axis
    name = TeX("Number of Crown-of-Thorns Starfish per Manta-Tow"))
ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_COTSwithdistlvl.png"), 
       plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
##########################################

##### reef_name single and compound #####
SingCompDistYears <- ifelse(is.na(obs_by_reef$single_or_compound), NA, obs_by_reef$YEAR)
ggplot(data = obs_by_reef,
       mapping = aes(x = YEAR, 
                     y = COVER)) +
  theme_classic() +
  geom_rect(xmin = SingCompDistYears, 
            xmax = ifelse(grepl("Unknown", obs_by_reef$recovYear), 
                          Inf,
                          as.numeric(obs_by_reef$recovYear)), 
            ymin = -Inf, ymax = Inf,
            mapping = aes(fill = obs_by_reef$single_or_compound),
            alpha = .3,
            na.rm = TRUE) + 
  scale_fill_manual(name = "Event Type",
                    values = c("orange", "yellow"),
                    na.translate = F) +
  geom_point(colour = "azure4",
             size = 1.5,
             alpha = 0.75) +
  geom_line(colour = "azure4",
            linewidth = 1,
            alpha = 0.75) +
  geom_vline(mapping = aes(xintercept = obs_by_reef$distYears), 
             colour = "gray40",
             na.rm = TRUE,
             linewidth = 1) +
  labs(x = "Year",
       y = "Cover (%)",
       colour = "Disturbance Type") + 
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015)) +
  scale_y_continuous(limits = c(0,50), 
                     breaks = seq(0, 50, 10))
if (isTimeBased) {
  ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_CompoundandSingleTimeBased.png"), 
         plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
} else {
  ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/", reef_name, "_CompoundandSingleRecovBased.png"), 
         plot = last_plot(), width=landscape_dims[1], height=landscape_dims[2])
}
##########################################

#### reef_name location (Crop after) #####
reef_coords_xy <- obs_by_reef$geometry %>% 
  st_coordinates()
leaflet() %>%
  setView(lng = mean(reef_coords_xy[,1]) - 1, 
          lat = mean(reef_coords_xy[,2]) + 1, 
          zoom = 5) %>%
  addProviderTiles("Esri.WorldGrayCanvas") %>%
  addPolygons(data = map_data, 
              popup = ~AREA_DESCR,
              color = "grey",
              opacity = 1,
              weight = 1) %>% 
  addCircleMarkers(lng = mean(reef_coords_xy[,1]),
                   lat = mean(reef_coords_xy[,2]),
                   color = "lightcoral",
                   opacity = 1,
                   weight = 4,
                   radius = 5,
                   fill = FALSE) %>%
  addLabelOnlyMarkers(lng = mean(reef_coords_xy[,1] + 1),
                      lat = mean(reef_coords_xy[,2]),
                      label = reef_name,
                      labelOptions = labelOptions(noHide = T, 
                                                  textOnly = TRUE, 
                                                  direction = "right",
                                                  style = list(
                                                    "color" = "lightcoral",
                                                    "font-weight" = "bold",
                                                    "font-size" = "18px"))) %>%
  addLabelOnlyMarkers(lng = 140,
                      lat = -25,
                      label = "Queensland",
                      labelOptions = labelOptions(noHide = T, 
                                                  textOnly = TRUE,
                                                  style = list(
                                                    "color" = "grey",
                                                    "font-weight" = "bold",
                                                    "font-size" = "12px"))) %>%
  saveWidget(paste0(reef_name, ".html"))
webshot(paste0(reef_name, ".html"), 
        paste0("../MPhil Thesis/Figures/Data Chapter/GBR_", reef_name, ".png"),
        vwidth = 650, vheight = 500, zoom = 8)
##########################################

##### VIS GRAVEYARD #####

# Geom point vis over map for recovery rate/time distribution, colour by number
# ggplot() +
#   theme_classic() +
#   geom_sf(data = map_data,
#           size = 0.5,
#           color = "azure4",
#           fill = NA) +
#   geom_point(data = all_reefs_sf[!is.na(all_reefs_sf$recovTime),],
#              aes(X, Y, colour = recovTime),
#              size = 1.5) +
#   labs(x = "Longitude",
#        y = "Latitude",
#        colour = "Average \nRecovery \nTime (years)")
# if (isTimeBased) {
#   ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/GBRAveRecoveryTime_TimeBased", 
#                 recov_yrs, "yr", mgmt_benefit, "mgmt.png"),
#          plot = last_plot(), width=5, height=5)
# } else {
#   ggsave(paste0("../MPhil Thesis/Figures/Data Chapter/GBRAveRecoveryTime_RecovBased", 
#                 recov_th, "th", mgmt_benefit, "mgmt.png"),
#          plot = last_plot(), width=5, height=5)
# }

## Plot the sectors with different colours
# ggplot() +
#   geom_sf(data = map_data,
#           mapping = aes(fill = AREA_DESCR), lwd = 0.01) +
#   theme_classic() +
#   labs(x = "Longitude",
#        y = "Latitude",
#        title = "The Great Barrier Reef",
#        subtitle = "by GBRMPA Management Area") +
#   scale_fill_brewer(name = "Region", palette = "Spectral") +
#   geom_sf(data = all_reefs_sf$geometry,
#           pch = 1,
#           color = "black")
# 
# # Leaflet plot of sectors
# leaflet() %>%
#   addTiles() %>%
#   addPolygons(data = map_data, 
#               popup = ~AREA_DESCR,
#               color = "blue",
#               opacity = 1,
#               weight = 1)
# 
# leaflet() %>%
#   addProviderTiles("Esri.WorldGrayCanvas") %>%
#   addPolygons(data = map_data, 
#               popup = ~AREA_DESCR,
#               color = "darkgrey",
#               opacity = 1,
#               weight = 1)
