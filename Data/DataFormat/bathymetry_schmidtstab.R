library(terra)
library(tidyterra)
library(sp)
library(tidyverse)

# Call in bathy data made by Sean ####

#If using the tif DEM -- too big to keep in GitHub
boysen_bathy <- rast("C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/BoysenRaster.tif")

#To plot contours easily all you have to do is:
contour(boysen_bathy) #this only works with the raster not the vector

depth_raster <- 1463.04 - boysen_bathy

# Replace all 0 values in depth raster with NA
depth_raster[depth_raster < 0.25] <- NA


# ggplot() +
#   geom_spatraster(data=depth_raster,aes(fill = lyr1)) +
#   scale_fill_viridis_c('Depth (m)', na.value = "transparent") +
#   theme_minimal() + 
#   guides(fill= guide_colorbar(reverse=T)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# ggsave('Figures/ASLO24/bathymetry.png',height=6.5,width=4.5,units='in',dpi=1200)
# 
# #dark theme 
# ggplot() +
#   geom_spatraster(data=depth_raster,aes(fill = lyr1)) +
#   scale_fill_viridis_c('Depth (m)', na.value = "transparent") +
#   dark_theme_minimal() + 
#   guides(fill= guide_colorbar(reverse=T)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# ggsave('Figures/ASLO24/bathymetry_dark.png',height=6.5,width=4.5,units='in',dpi=1200)


#If using the contour lines shapefile
# boysen_bathy <- vect("C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/Boysen_Contours/Boysen_Contours/All_Boysen_Contours.shp")
#You can plot either of these easily using base plot or use the tidyterra package for ggplot


#As far as cross sectional analysis you'll probably want to use the raster from the DEM tif. Terra has awesome documentation that gives a great idea about the tools available to you, but if you want to talk more in person about what you want to do I'd be happy to.


# from GPT

# Convert elevation to positive depth values
depth <- abs(boysen_bathy)

# Calculate cross-sectional areas
areas <- cellSize(depth)

# Extract the depth and area values from the SpatRaster object
depth_values <- values(depth)
area_values <- values(areas)

# Combine the depth and area values into a dataframe
depth_area_df <- data.frame(
  Depth = depth_values,
  Area = area_values
) |> distinct() 

head(depth_area_df)

depth_area_aggregated <- depth_area_df |>
  mutate(depth_m = 1463.04 - Depth.lyr1) |> # elevation of reservoir minus elevation of DEM
  group_by(depth_m) |>
  summarise(area_m2 = sum(area)) |>
  ungroup() |>
  mutate(depth_m_rounded = plyr::round_any(depth_m,0.5)) |>
  group_by(depth_m_rounded) |>
  summarise(area_m2_rounded = sum(area_m2)) |>
  ungroup()

hypso <- depth_area_aggregated |>
  rename(depths=depth_m_rounded,
         areas=area_m2_rounded)



# # Load your points or polygons (replace "points.shp" with your file)
# samplingsites <- BoysenChem |>
#   select(WaterbodyName, Latitude, Longitude) |>
#   distinct() 
# coordinates(samplingsites) <- c("Longitude", "Latitude")
# 
# points <- vect(samplingsites)
# 
# # Extract depth measurements at points
# depth_measurements <- terra::extract(depth, points)
# areas_at_points <- terra::extract(areas, points)





# Find Schmidt's Stability ####
source('Data/CALL_DATA_LIB.R')
library(rLakeAnalyzer)

hypso <- read.csv('Data/Simplified_bathymetry.csv')


smidts_stability <- BoysenProfile |>
  drop_na() |>
  mutate(depth = plyr::round_any(depth_m,0.5)) |>
  group_by(WaterbodyName, depth, CollDate) |>
  summarise(temp=mean(temp_C))|>
  # left_join(hypso, by=c('depth_m'= 'depths')) |>
  # drop_na(areas) |>
  group_by(WaterbodyName, CollDate) |>
  summarise(SS = schmidt.stability(temp, depth, hypso$areas, hypso$depths, sal=0)) |>
  as.data.frame() |>
  mutate(SS = ifelse(SS<0,0,SS)) |> # two points where stability is <0, that defies laws of physics so make 0
  mutate(Year = year(CollDate),
         Month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  group_by(Year,Month) |>
  mutate(average_mon_whole_SS = mean(SS)) 

ggplot(smidts_stability) +
  geom_point(aes(Month, SS, color=WaterbodyName)) +
  geom_point(aes(Month, average_mon_whole_SS),size=4,shape=21) +
  scale_color_viridis_d('', option='turbo') +
  theme_bw() +
  facet_wrap(~Year) +
  labs(x='',y='Schmidt Stability Index'~(J~m^-2))


ggplot(smidts_stability) +
  geom_boxplot(aes(WaterbodyName, SS, group=WaterbodyName))+
  labs(x='',y='Schmidt Stability Index'~(J~m^-2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))



write.csv(smidts_stability, 'Data/Schmidts_Stability.csv')
write.csv(hypso, 'Data/Simplified_bathymetry.csv')



# try making a 3d map
# library(marmap)
# bathy_xyz <- terra::as.data.frame(boysen_bathy,xy=TRUE)
# head(bathy_xyz)
# 
# bath <- marmap::as.bathy(bathy_xyz |> dplyr::select(1:3))
# 
# library(lattice)
# lattice::wireframe(unclass(bath), shade = TRUE, aspect = c(1/2, 0.1))

# try making a cross section
# define transect
start_point <- c(lon = -108.180325, lat = 43.173358) 
end_point <- c(lon = -108.180325, lat = 43.416572) 

#extract data along transect
# Create a sequence of points along the transect
lon_seq <- seq(start_point["lon"], end_point["lon"], length.out = 5000)
lat_seq <- seq(start_point["lat"], end_point["lat"], length.out = 5000)


# Create a data frame of these points
transect_points <- data.frame(lon = lon_seq, lat = lat_seq)

# Extract depth values at these points
transect_points$depth <- terra::extract(boysen_bathy, cbind(transect_points$lon, transect_points$lat))[,1]



# Create the cross-section plot
ggplot(transect_points, aes(x = lat, y = 1463.04 -depth)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  theme_minimal() +
  labs(x = "Latitude", y = "Depth (m)", title = "Lake Bathymetry Cross-Section") +
  scale_y_reverse()
