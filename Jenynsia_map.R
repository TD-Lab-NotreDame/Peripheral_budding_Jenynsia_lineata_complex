#libraries
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)

# code to make whole of South America with box around study area
south_america_map <- map_data("world", region = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador", 
                                                  "Guyana", "Paraguay", "Peru", "Suriname", "Uruguay", "Venezuela", 
                                                  "French Guiana"))

# Define your bounding box coordinates
longitude_min <- -78
longitude_max <- -40
latitude_min <- -15
latitude_max <- -40

# Create a data frame for the polygon
bounding_box <- data.frame(
  longitude = c(longitude_min, longitude_max, longitude_max, longitude_min, longitude_min),
  latitude = c(latitude_min, latitude_min, latitude_max, latitude_max, latitude_min)
)


ggplot() +
  geom_polygon(data = south_america_map, aes(x = long, y = lat, group = group), 
               fill = "lightgrey", color = "black") +  # Draw a red dashed line at -40° S
  coord_map("mercator") +  # Apply the Mercator projection
  xlim(-88, NA) +
  geom_polygon(data = bounding_box, aes(x = longitude, y = latitude), fill = NA, color = "skyblue4", size = 1) +  # Adjust limits for better visualization
  theme_void()
##########
#code to make zoomed in map around study area and add sampling points
#sampling locations
data<-read_csv("data/jenynsia_data_mapping.csv")
#map data
world <- ne_countries(scale="medium", returnclass = "sf")

# bounding box of interest
bbox <- c(xmin = -80, xmax = -40, ymin = -40, ymax = -10)
#get river information 
rivers <- ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")
# Crop to bounding box of your South America map
rivers_cropped <- st_intersection(rivers, st_bbox(south_america_filtered) %>% st_as_sfc())

# Crop the world map data to this bounding box
south_america_filtered <- st_crop(world, bbox)

#plot
ggplot(data = south_america_filtered) +
  geom_sf(fill = "lightgrey", color = "black", size = 0.3) +  # Custom fill and border color
  geom_sf(data = rivers_cropped, color = "lightblue3", size = 0.3)+# adds rivers
  coord_sf(expand = FALSE) +  # Prevents padding around the map
  theme_minimal() +  # Minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  xlab("Longitude") +  geom_point(data = data, aes(x = long, y = lat, fill = population_map, shape= population_map), size = 2) +
  scale_fill_manual(
    values = c("#EE9A00", "#6CC7AB", "#B95251", "#4D4D4D"))+
  scale_shape_manual(
    values = c(22, 21, 23, 24))+
  labs(fill = NULL) +
  labs(shape = NULL) +
  scale_x_continuous(breaks = seq(-80, 40, by = 10)) +
  ylim(-40, -15) +
  xlab("Longitude") + ylab("Lattitude") +
  theme_minimal() +theme(legend.position =  c(0.4, -.25),  legend.direction = "horizontal") +
  ylab("Latitude") +
  annotation_scale(location = "bl", width_hint = 0.2) +  # Adds scale bar
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +  # Adds north arrow
  theme(# Removes grid lines
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")  # Adds margin around the plot
  )

