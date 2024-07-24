#--------------------------------------------#
# Nutrient concentration maps
#--------------------------------------------#

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')
library(sf)
shapefile <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')


# 2. Look at nutrients by site/month ####
ggplot(BoysenNutrient |> filter(ShortName_Revised=='TN')) +
  geom_jitter(aes(month, ChemValue,color=WaterbodyName)) +
  facet_wrap(~Year)

ggplot(BoysenNutrient |> filter(ShortName_Revised=='TP')) +
  geom_jitter(aes(month, ChemValue,color=WaterbodyName)) +
  facet_wrap(~Year)


# 3. Calculate average nutrients by site ####
ave_nutrient <- BoysenNutrient |>
  filter(ShortName_Revised%in%c('TN','TP')) |>
  group_by(WaterbodyName, ShortName_Revised, Longitude, Latitude) |>
  summarise(meanconc=mean(ChemValue)) |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 

ave_cyano <- BoysenPhyto_cat |>
  left_join(BoysenNutrient |> select(WaterbodyName,Longitude,Latitude) |> distinct()) |>
  group_by(WaterbodyName, Longitude, Latitude) |>
  summarise(meanconc=mean(Cyanobacteria))  |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 


# 4. Plot average nutrients by site ####
ggplot()+
  geom_sf(shapefile, mapping=aes(), fill='black',color='white') +
  geom_sf(ave_nutrient |> filter(ShortName_Revised=='TN'),mapping=aes(size=meanconc),fill='red4',shape=21) +
  dark_theme_void() +
  labs(x='',y='') +
  theme(legend.title = element_blank())
ggsave('Figures/ASLO24/TN_map.png',height=6.5,width=4.5, units='in',dpi=1200)

ggplot()+
  geom_sf(shapefile, mapping=aes(), fill='black',color='white') +
  geom_sf(ave_nutrient |> filter(ShortName_Revised=='TP'),mapping=aes(size=meanconc),fill='#336a98',shape=21) +
  dark_theme_void() +
  labs(x='',y='') +
  theme(legend.title = element_blank())
ggsave('Figures/ASLO24/TP_map.png',height=6.5,width=4.5, units='in',dpi=1200)


ggplot()+
  geom_sf(shapefile, mapping=aes(), fill='black',color='white') +
  geom_sf(ave_cyano,mapping=aes(size=meanconc)) +
  dark_theme_void() +
  labs(x='',y='') +
  theme(legend.title = element_blank())
ggsave('Figures/ASLO24/cyano%_map.png',height=6.5,width=4.5, units='in',dpi=1200)
