library(terra)
library(tidyterra)
library(sp)



#If using the tif DEM
boysen_bathy <- rast("C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/BoysenRaster.tif")

#To plot contours easily all you have to do is:
contour(boysen_bathy) #this only works with the raster not the vector

#If using the contour lines shapefile
#boysen_bathy <- vect("C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/Boysen_Contours/Boysen_Contours/All_Boysen_Contours.shp")
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
  summarise(area_m2_rounded = sum(area_m2))

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





# attempt some things using hypsometry data now
source('Data/CALL_DATA_LIB.R')
library(rLakeAnalyzer)


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



