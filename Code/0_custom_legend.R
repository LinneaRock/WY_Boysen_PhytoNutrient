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


# custom legend
ggplot()+
  geom_sf(shapefile, mapping=aes(),fill='white',color='grey20') +
  geom_sf(sites_sf,mapping=aes(fill=WaterbodyName),size=5, shape=21) +
  # geom_text(sites_sf, mapping=aes(lon,lat, label=WaterbodyName),nudge_x=0,hjust=1.05) +
  geom_text(sites_sf, mapping=aes(lon,lat, label=WaterbodyName),vjust=2,size=5) +
  theme_void() +
  labs(x='',y='') +
  scale_fill_viridis_d('',option='magma') +
  theme(legend.position='none') +
  expand_limits(x=-108)

# dark theme 
ggplot()+
  geom_sf(shapefile, mapping=aes(), fill='black',color='white') +
  geom_sf(sites_sf,mapping=aes(fill=WaterbodyName),size=5,shape=21) +
  # geom_text(sites_sf, mapping=aes(lon,lat, label=WaterbodyName),nudge_x=0,hjust=1.05) +
  geom_text(sites_sf, mapping=aes(lon,lat, label=WaterbodyName),vjust=2,size=5) +
  dark_theme_void() +
  labs(x='',y='') +
  scale_fill_viridis_d('',option='magma') +
  theme(legend.position='none') +
  expand_limits(x=-108)


# tribs custom legend
tribs_sf <- read.csv('Data/trib_load_flux.csv') |>
  dplyr::select(WaterbodyName, Longitude, Latitude) |>
  distinct() |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 

tribs <- read.csv('Data/trib_load_flux.csv') |>
  dplyr::select(WaterbodyName, Longitude, Latitude) |>
  distinct() 

trib_colors <- c('#332288','#882255','#E6AA68','#DDCC77')

# flowlines
NHDflowline <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/NHD/NHDPLUS_H_1008_HU4_GDB.gdb', layer = 'NHDFlowline')
st_crs(NHDflowline)
NHDflowline <- st_transform(NHDflowline, crs=4326)

st_crs(NHDflowline)==st_crs(tribs_sf)

tribs_flowline <- NHDflowline |>
  filter(GNIS_Name %in% c('Fivemile Creek', 'Muddy Creek', 'Wind River'))

# Drop the Z and M dimensions from the geometries
tribs_flowline <- st_zm(tribs_flowline, drop = TRUE)

tribs_sf$WaterbodyName <- factor(tribs_sf$WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet'))

ggplot()+
  geom_sf(tribs_flowline, mapping=aes()) +
  geom_sf(shapefile, mapping=aes(),fill='white',color='grey20') +
  geom_sf(tribs_sf, mapping=aes(fill=WaterbodyName),size=5,shape=23,color='white') +
  geom_text(tribs, mapping=aes(Longitude, Latitude, label=WaterbodyName),vjust=-1, hjust=1,size=5) +
  scale_fill_manual('',values =trib_colors) +
  theme_void() +
  labs(x='',y='') +
  theme(legend.position = 'none') +
  coord_sf(xlim = c(-108.5, -108), ylim = c(43.1, 43.5)) 

# library(export)
# graph2ppt(file='Figures/custom_legend_tribs.pptx', append=TRUE)


