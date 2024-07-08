############################
# Figure 1 a beautiful map 
############################

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')
library(sf)
library(maps)
library(ggspatial)
library(raster)

# lake shapefile 
lake_shapefile <- st_read('C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')
st_crs(lake_shapefile)
lake_shapefile <- st_transform(lake_shapefile, crs = 4326)

# lake sites
sites_sf <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                  fileEncoding="latin1") |> # doesn't matter which csv we call
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate)) |>
  filter(Year>=2020) |>
  dplyr::select(WaterbodyName, Latitude, Longitude) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  distinct() |>
  mutate(WaterbodyName=trimws(WaterbodyName),
         lat=Latitude,
         lon=Longitude) |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 
st_crs(sites_sf)

# WY state outline
WY <- map_data('state', region='Wyoming')


# DEM 
dem_raster <- raster('C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/WY_elevation')
# Create a black-and-white color palette with darker colors for higher elevations
bw_palette <- gray.colors(100, start = 0, end = 1, gamma = 2.2, rev = TRUE)
plot(dem_raster, col = bw_palette, main = "Digital Elevation Model with Custom Colors")
# thx gpt!! XD
st_crs(dem_raster)
dem_raster <- projectRaster(dem_raster, crs=4326)

# watershed shapefile
ws_shp <- st_read('Data/Watershed_shapefile/area-of-interest.shp')
ws_shp <- st_transform(ws_shp, crs=4326)


# flowlines
NHDflowline <- st_read('C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/NHD/NHDPLUS_H_1008_HU4_GDB.gdb', layer = 'NHDFlowline')
st_crs(NHDflowline)
NHDflowline <- st_transform(NHDflowline, crs=4326)

# tributary/outlet locations
tribs_sf <- BoysenTribs |>
  dplyr::select(WaterbodyName, Latitude, Longitude) |>
  distinct() |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 
st_crs(tribs_sf)


#bathymetry 
boysen_bathy <- st_read("C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/Boysen Shapefile/boysen_Shapefile.shp")
st_crs(boysen_bathy)

# 2. check projections ####
st_crs(lake_shapefile)==st_crs(sites_sf)
st_crs(sites_sf)==st_crs(dem_raster)
st_crs(dem_raster)==st_crs(ws_shp)
st_crs(ws_shp)==st_crs(NHDflowline)
st_crs(NHDflowline)==st_crs(tribs_sf)
st_crs(tribs_sf)==st_crs(boysen_bathy)


# 3. Crop any spatial opjects as necessary ####
# cropping DEM 
cropped_dem <- terra::crop(dem_raster, ws_shp, mask=T)
str(cropped_dem)
st_crs(cropped_dem)


# cropping flowlines
cropped_flowline <- st_crop(NHDflowline, ws_shp)
croppe <- st_crop(NHDflowline, xmin=-108.339946,ymin=43.12895,xmax=-108.039283,ymax=43.427151)

library(mapview)
mapview(croppe)

# need to do this in gis


main_IDs <- c(23002600043697)






plot(dem_raster)
#plot(cropped_flowline, col='grey',add=T, alpha=0.1)
plot(ws_shp, col='transparent', add=T)
plot(lake_shapefile, col='black', add=T)
plot(sites_sf, col='red', add=T)


plot(cropped_dem) 
plot(ws_shp)
plot(NHDflowline)

ggplot() +
  geom_sf(NHDflowline, mapping=aes())
