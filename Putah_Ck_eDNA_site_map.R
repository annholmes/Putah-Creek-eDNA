#Putah Creek eDNA study (Holmes et al. 2025)

#R script for a map of lower Putah Creek with a call out map of California
#eDNA sample sites sites sampled 2019-2022 are labeled

library(sf) #import shapefiles
library(ggplot2) 
library(tidyverse)
library(leaflet) #mapping
library(htmlwidgets) #workaround to save an interactive leaflet map
library(webshot) #saves a static image from leaflet map
#library(ggimage) #combine maps by inserting leaflet map as an image over regional map
library(magick)
library(grid)
library(cowplot) #arrange multiple plots

setwd("~/PC") #set user specific working directory
options(scipen = 999) #so numbers don't display as scientific notation

#shapefile for map of California
#downloaded from Stanford Geospatial Center library 
#https://library.stanford.edu/libraries/stanford-geospatial-center
##all other files must remain in the folder to read the shape file
NAmer <- st_read("~/PC/stanford-ns372xw1938-shapefile/ns372xw1938.shp", 
                 stringsAsFactors = FALSE) 

#create an 
cali <-ggplot() +
  geom_sf(data = NAmer, fill = "#ECFADC") +
  theme_bw()+
  coord_sf(default_crs = sf::st_crs(4326),
           xlim = c(-125, -115), 
           ylim = c(32, 44)) +
  scale_x_discrete(breaks = c(-125, -120, -115), 
                   labels = c("-125", "-120", "-115")) + 
  scale_y_discrete(breaks = c(32, 36, 40, 44), 
                   labels = c("32", "36", "40", "44")) +
  annotate("rect", xmin=c(-122), xmax=c(-121.5), 
           ymin=c(38.4) , ymax=c(38.6), 
           alpha=0.2, color="blue")+
  annotate(geom = "text", x = -122.5, y = 32.5, label = "Pacific Ocean", 
           fontface = "italic", color = "grey22", size = 3)+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())

cali

#full list of 11 eDNA sample sites
points_df = tribble(~name, ~lat, ~lng,
                    'SCCA', 38.5128,-122.0972,
                    'WPCP', 38.5204,-121.9670,
                    'I505', 38.5261,-121.9513,
                    'SW1D', 38.5325,-121.9355,
                    'SW2D', 38.5313,-121.9313,
                    'RRAN', 38.5379, -121.8650,
                    'PEDR', 38.5268,-121.8035,
                    'APOP', 38.5242,-121.7889,
                    'OLDR', 38.5172,-121.7562,
                    'MACU', 38.5176,-121.6967,
                    'MACD', 38.5192, -121.6967)
# Convert the data frame to a spatial points data frame 
Sites = st_as_sf(points_df,  
                 coords = c("lng", "lat"), crs = 4326) 

#create custom data frames for adjustments to locations of site labels relative to points
#adjusted so labels of sites that are close by don't overlap
points_df_top = tribble(~name_top, ~lat, ~lng,
                        'SCCA', 38.5128,-122.0972,
                        'SW1D', 38.5325,-121.9355,
                        'RRAN', 38.5379, -121.8650,
                        'PEDR', 38.5268,-121.8035,
                        'OLDR', 38.5172,-121.7562)

points_df_I505 = tribble(~name_I505, ~lat, ~lng,
                         'I505', 38.5261,-121.9513)

points_df_bottom = tribble(~name_bottom, ~lat, ~lng,
                           'WPCP', 38.5204,-121.9670)

points_df_APOP = tribble(~name_APOP, ~lat, ~lng,
                         'APOP', 38.5242,-121.7889)

points_df_left = tribble(~name_left, ~lat, ~lng,
                         'MACU', 38.5176,-121.6967)

points_df_right = tribble(~name_right, ~lat, ~lng,
                          'SW2D', 38.5313,-121.9313,
                          'MACD', 38.5192, -121.6930)

pc <-leaflet(Sites) %>%  
  addProviderTiles('Esri.WorldTopoMap') %>% 
  addCircles(points_df_top$lng, points_df_top$lat, 
             color = "black",
             opacity = 1,
             label = points_df_top$name_top, 
             labelOptions = labelOptions(noHide = T, 
                                         direction = "top",
                                         offset = c(0, 8),
                                         textOnly = TRUE,
                                         textsize = "15px"),
             group = "Show All")%>%
  addCircles(points_df_bottom$lng, points_df_bottom$lat, 
             color = 'black',
             opacity = 1,
             label = points_df_bottom$name_bottom, 
             labelOptions = labelOptions(noHide = T, 
                                         direction = "bottom",
                                         offset = c(0, -8),
                                         textOnly = TRUE,
                                         textsize = "15px"),
             group = "Show All")%>%
  addCircles(points_df_right$lng, points_df_right$lat, 
             color = 'black',
             opacity = 1,
             label = points_df_right$name_right, 
             labelOptions = labelOptions(noHide = T, 
                                         direction = "right",
                                         offset = c(6, -6),
                                         textOnly = TRUE,
                                         textsize = "15px"),
             group = "Show All")%>%
  addCircles(points_df_left$lng, points_df_left$lat, 
             color = 'black',
             opacity = 1,
             label = points_df_left$name_left, 
             labelOptions = labelOptions(noHide = T, 
                                         direction = "left",
                                         offset = c(-6, 8),
                                         textOnly = TRUE,
                                         textsize = "15px"),
             group = "Show All")%>%
  addCircles(points_df_I505$lng, points_df_I505$lat, 
             color = 'black',
             opacity = 1,
             label = points_df_I505$name_I505, 
             labelOptions = labelOptions(noHide = T, 
                                         direction = "top",
                                         offset = c(-15, 10),
                                         textOnly = TRUE,
                                         textsize = "15px"),
             group = "Show All")%>%
  addCircles(points_df_APOP$lng, points_df_APOP$lat, 
             color = 'black',
             opacity = 1,
             label = points_df_APOP$name_APOP, 
             labelOptions = labelOptions(noHide = T, 
                                         direction = "bottom",
                                         offset = c(-10, -10),
                                         textOnly = TRUE,
                                         textsize = "15px"),
             group = "Show All")%>%
  addScaleBar(position = "bottomright")

## save html to png
saveWidget(pc, "temp.html", selfcontained = FALSE)
pc_map <-webshot("temp.html", 
                 file = "Putah_Ck_map.png",
                 cliprect = "viewport",
                 zoom = 4) #zoom in for better image quality

#create object from saved image
pc_map_crop <-image_read("~/PC/Putah_Ck_map_cropped.png")

#create more room on regional map for adding the Putah Creek map
cali_roomy <-plot_grid(cali,"","",
                       ncol = 3)

#combine regional California map with Putah Creek eDNA sampling site map
final_map <- ggdraw() +
  draw_plot(cali_roomy) +
  draw_image(pc_map_crop,  x = 0.162, y = 0.02, 
             scale = .9) 

#save plot with call out map on left and site map of Putah Creek on right
save_plot("Putah_Ck_eDNA_final_map.png", 
          final_map, 
          dpi = 300,
          base_height = 6, #height in inches, equal to 15.24cm
          base_width = 12, #height in inches, equal to 30.48cm
          bg = "white")
