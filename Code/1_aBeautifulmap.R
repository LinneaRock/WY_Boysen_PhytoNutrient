#~~~~~~~~~~~~~~~~#
# A beautiful map 
#~~~~~~~~~~~~~~~~#

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')
library(sf)
#library(maps)
library(ggspatial)
library(raster)
library(usmap)

# lake shapefile 
lake_shapefile <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')
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

# # WY state outline
# WY <- map_data('state', region='Wyoming')
# # U.S. outline 
# us <- map_data('usa') #|>
#   #st_as_sf(coords=c('long','lat'), crs=4326) 


# DEM 
dem_raster <- raster('C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/WY_elevation')
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
NHDflowline <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/NHD/NHDPLUS_H_1008_HU4_GDB.gdb', layer = 'NHDFlowline')
st_crs(NHDflowline)
NHDflowline <- st_transform(NHDflowline, crs=4326)

# tributary/outlet locations
tribs_sf <- BoysenTribs |>
  dplyr::select(WaterbodyName, Latitude, Longitude) |>
  distinct() |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 
st_crs(tribs_sf)

tribs_sf$WaterbodyName <- factor(tribs_sf$WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet'))

#bathymetry 
# boysen_bathy <- st_read("C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/Boysen Shapefile/boysen_Shapefile.shp")
# st_crs(boysen_bathy)

# 2. check projections ####
st_crs(lake_shapefile)==st_crs(sites_sf)
st_crs(sites_sf)==st_crs(dem_raster)
st_crs(dem_raster)==st_crs(ws_shp)
st_crs(ws_shp)==st_crs(NHDflowline)
st_crs(NHDflowline)==st_crs(tribs_sf)
#st_crs(tribs_sf)==st_crs(boysen_bathy)


# 3. Crop any spatial opjects as necessary ####
# cropping DEM 
cropped_dem <- terra::crop(dem_raster, ws_shp, mask=T)
str(cropped_dem)
st_crs(cropped_dem)


# cropping flowlines
cropped_flowline <- st_intersection(NHDflowline, ws_shp)

cropped_flowline_majors <- cropped_flowline |>
  drop_na(GNIS_Name)


# 4. check out objects with mapview ####
library(mapview)
mapview(ws_shp) + mapview(cropped_flowline_majors,color='black') + mapview(tribs_sf,color='red') + mapview(lake_shapefile, color='black')

mapview(dem_raster) + mapview(ws_shp) + mapview(cropped_flowline_majors,color='black') + mapview(tribs_sf,color='red') + mapview(lake_shapefile, color='black')

mapview(cropped_dem) + mapview(ws_shp) + mapview(cropped_flowline_majors,color='black') + mapview(tribs_sf,color='red') + mapview(lake_shapefile, color='black')


# 5. Part A - Watershed ####
cropped_dem_spdf <- as(cropped_dem, "SpatialPixelsDataFrame")
cropped_dem_df <- as.data.frame(cropped_dem_spdf)
colnames(cropped_dem_df) <- c("value", "x", "y")

trib_colors <- c('#332288','#882255','#E6AA68','#DDCC77')

library(ggnewscale)
a<-ggplot() +
  geom_tile(cropped_dem_df, mapping=aes(x,y, fill=value), alpha=0.85) +
  scale_fill_gradientn('Elevation (m)', colors=bw_palette) +
  #scale_fill_viridis_c('Elevation (m)', option='rocket', direction=-1) +
  geom_sf(ws_shp, mapping=aes(), color='black', fill=NA, linewidth=0.75) +
  geom_sf(cropped_flowline_majors, mapping=aes(), color='black',linewidth=0.25) +
  new_scale_fill() +
  geom_sf(tribs_sf, mapping=aes(fill=WaterbodyName),size=2,shape=23,color='black') +
  scale_fill_manual('',values =trib_colors) +
  geom_sf(lake_shapefile, mapping=aes(), fill='black') +
  theme_minimal() +
  labs(x='',y='')

# 6. Part A, inset ####
ws_shp_transformed <- usmap_transform(ws_shp)

ainsert <- plot_usmap(exclude = c("AK","HI"), color='gray50') +
  geom_sf(ws_shp_transformed, mapping=aes(), color='black', fill='black')


# 7. Part B - Lake ####
b <- ggplot()+
  geom_sf(lake_shapefile, mapping=aes(),fill='white',color='grey20') +
  geom_sf(sites_sf,mapping=aes(fill=WaterbodyName),size=4, shape=21) +
  theme_minimal() +
  labs(x='',y='') +
  scale_fill_viridis_d('',option='magma') +
  theme(axis.text.x=element_text(angle=45))

# 8. Save plots and manually manipulate in PP :( ####  
a 
ggsave('Figures/watershed.png', width=4.5,height=6.5,units='in',dpi=1200)
a + theme(legend.position = 'none')
ggsave('Figures/watershed_nolegend.png', width=4.5,height=6.5,units='in',dpi=1200)

ainsert
ggsave('Figures/countryinsert.png', width=4.5,height=6.5,units='in',dpi=1200)

b
ggsave('Figures/boysen.png', width=4.5,height=6.5,units='in',dpi=1200)
