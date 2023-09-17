# Packages
library(hexbin)
library(RColorBrewer)

compound_vals <- st_transform(all_samples, crs = st_crs(sector_boundaries))

x <- st_coordinates(compound_vals$geometry)[,1]
y <- st_coordinates(compound_vals$geometry)[,2]

compound_plot <- ggplot() +
  geom_sf(data = sector_boundaries, lwd = 0.01) +
  theme_classic() +
  labs(x = "Longitude",
       y = "Latitude",
       title = "Sample reefs chosen to manage in a compound disturbance system",
       subtitle = "in the Great Barrier Reef by GBRMPA Management Area") +
  geom_hex(compound_vals, 
           mapping = aes(x, y),
           binwidth = c(0.5,0.5)) + 
  labs(fill = "Count")
compound_plot
