#--------------------------------------------#
# Linnea's custom Boysen Reservoir map legend
#--------------------------------------------#

# 1. load libraries, data ####
library(sf)
library(tidyverse)

shapefile <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')

sites <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                   fileEncoding="latin1") |> # doesn't matter which csv we call
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate)) |>
  filter(Year>=2020) |>
  select(WaterbodyName, Latitude, Longitude) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  distinct()

sites_sf <- sites |>
  mutate(WaterbodyName=trimws(WaterbodyName),
         lat=Latitude,
         lon=Longitude) |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 



ggplot()+
  geom_sf(shapefile, mapping=aes(),fill='white',color='grey20') +
  geom_sf(sites_sf,mapping=aes(color=WaterbodyName),size=5) +
  # geom_text(sites_sf, mapping=aes(lon,lat, label=WaterbodyName),nudge_x=0,hjust=1.05) +
  geom_text(sites_sf, mapping=aes(lon,lat, label=WaterbodyName),vjust=2,size=5) +
  theme_void() +
  labs(x='',y='') +
  scale_color_viridis_d('',option='turbo') +
  theme(legend.position='none') +
  expand_limits(x=-108)



