#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Schmidt's Stability Figure 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Code to calculation Schmidt's Stability is located in Data/DataFormat/bathymetry_schmidtstab.R


library(tidyverse)
library(terra)
library(tidyterra)
library(sp)
library(sf)
library(patchwork)

# read in stability data 
SS <- read.csv('Data/Schmidts_Stability.csv') |>
  rename(month=Month)

# determine how different sampling locations are in terms of stability 
h <- aov(SS~month, SS)
tukey <- TukeyHSD(h)
library(multcompView)
cld <- multcompLetters4(h, tukey)
cld2 <- data.frame(letters = cld$month$Letters)
cld2$month <- rownames(cld2)
sig.letters <- cld2

# significance for plotting
sig.letters <- sig.letters |>
  drop_na(letters) 

means <- left_join(SS, sig.letters) |>
  group_by(month, letters) |>
  summarise(max.result = max(SS, na.rm = TRUE)) |>
  distinct()

# make SS plot
ssplot <- SS |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  mutate(month=factor(month, levels=c('May','Jun','Jul','Aug','Sep','Oct'))) |>
  ggplot() +
  geom_boxplot(aes(month, SS)) +
  geom_jitter(aes(month, SS, color=WaterbodyName),alpha=0.5) +
  scale_color_viridis_d('', option='turbo') +
  geom_text(means, mapping=aes(month, 
                               max.result+100, label = letters), 
            size=4) +
  labs(x='',y='Schmidt Stability Index'~(J~m^-2)) +
  theme_minimal() +
  theme(legend.position = 'none')



# add bathy plot 
# create sites sf dataframe
sites_sf <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                     fileEncoding="latin1") |> # doesn't matter which csv we call
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate)) |>
  filter(Year>=2020) |>
  select(WaterbodyName, Latitude, Longitude) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  distinct() |>
  mutate(WaterbodyName=trimws(WaterbodyName),
         lat=Latitude,
         lon=Longitude) |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 

# read in tif DEM -- too big to keep in GitHub
boysen_bathy <- rast("C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/BoysenRaster.tif")

depth_raster <- 1463.04 - boysen_bathy

# Replace all 0 values in depth raster with NA
depth_raster[depth_raster < 0.25] <- NA

# plot
bathy_plot <- ggplot() +
  geom_spatraster(data=depth_raster,aes(fill = lyr1)) +
  geom_sf(sites_sf, mapping=aes(color=WaterbodyName), size=3) +
  scale_fill_gradient('Depth (m)', na.value = "transparent", high='grey90', low='grey10') +
  scale_color_viridis_d('', option='turbo') +
  theme_void() +
  guides(fill= guide_colorbar(reverse=T)) 


bathy_plot + ssplot +
  plot_annotation(tag_levels = 'a',
                  tag_suffix = ')')
ggsave('Figures/bathy_stability.png',height=4.5, width=8.5, units='in',dpi=1200)

