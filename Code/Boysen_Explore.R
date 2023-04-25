#-----------------------------------------------#
# Exploratory plots for April 2023 meeting
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')

params <- unique(ChemPhys$ShortName_Revised) # which parameters are available in this dataset? 
params

# 1. Prep data for plotting ####
BoysenChemPhys <- ChemPhys |>
  select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  mutate(ChemUnits = sub('/l', '/L', ChemUnits)) |>
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  select(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, BelowDet, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct()

## how many of each parameter are below detection at each location? ####
n_belowDet <- BoysenChemPhys |>
  group_by(WaterbodyName, ShortName_Revised) |>
  add_count() |>
  mutate(range = paste(min(ChemValue), max(ChemValue), sep='-')) |>
  ungroup() |>
  select(WaterbodyName, ShortName_Revised, BelowDet, n, range) |>
  group_by(WaterbodyName, ShortName_Revised, n, range) |>
  summarise(belowDet = sum(BelowDet),
            percLow = belowDet/n *100) |>
  distinct()

# 2. Basic timeseries! ####
params <- unique(BoysenChemPhys$ShortName_Revised)

for(p in 1:length(params)) {

    ggplot(BoysenChemPhys |> filter(ShortName_Revised==params[p])) + 
      geom_point(aes(CollDate, SampleDepth, color=ChemValue)) +
      labs(y='Sample depth (m)', x='', title=params[p]) +
      scale_y_reverse() +
      facet_wrap(~WaterbodyName, scales='free_y') +
      scale_color_viridis_c(paste(params[p],unique((BoysenChemPhys |>
                                                      filter(ShortName_Revised==params[p]))$ChemUnits), 
                                  sep = ' ')) +
      theme_classic() 
    ggsave(paste0('Figures/Boysen_explore/',params[p],'_depth.png'),width=10, height=8, units='in')
  
    ggplot(BoysenChemPhys |> filter(ShortName_Revised==params[p],
                                    SampleDepth <= 3)) + 
      geom_point(aes(CollDate, ChemValue)) +
      labs(y=paste(params[p],unique((BoysenChemPhys |> 
                                       filter(ShortName_Revised==params[p]))$ChemUnits), 
                   sep = ' '), x='') +
      facet_wrap(~WaterbodyName, scales='free_y') +
      theme_classic() 
    ggsave(paste0('Figures/Boysen_explore/',params[p],'_epi.png'),width=10, height=8, units='in')
    
}

# 3. Nutrient timeseries including all nuts ####
# epilimnion only
ggplot(BoysenChemPhys |> filter(SampleDepth <= 3,
                                ShortName_Revised != 'Chlorophyll a (phytoplankton)')) + 
  geom_point(aes(CollDate, ChemValue, color=ShortName_Revised)) +
  facet_wrap(~WaterbodyName, scales='free_y') +
  theme_classic() 

# 4. What frequency are the data collected at each location? ####
boysensites <- unique(BoysenChemPhys$WaterbodyName) # total of 13 sites 

Annualfreq <- BoysenChemPhys |>
  filter(ShortName_Revised != 'Total Kjeldahl Nitrogen (unfiltered)') |>
  group_by(WaterbodyName, Latitude, Longitude, year(CollDate)) |>
  add_count() |>
  ungroup() |>
  select(-StationID, -BelowDet, -ChemValue, -ChemUnits, -SampleDepth, -CollDate,-ShortName_Revised) |>
  unique()

SiteVisitsMonthly <- BoysenChemPhys |>
  filter(ShortName_Revised != 'Total Kjeldahl Nitrogen (unfiltered)') |>
  select(-StationID, -Latitude, -Longitude, -BelowDet, -ChemValue, -ChemUnits, -SampleDepth) |>
  mutate(year = year(CollDate),
         month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  select(-CollDate) |>
  distinct() |>
  group_by(year, month, ShortName_Revised) |>
  add_count() |>
  ungroup() |>
  distinct() |>
  select(-ShortName_Revised, -WaterbodyName) |>
  distinct()

### many sites per monthly visit over time ####
ggplot(SiteVisitsMonthly, aes(month, year, fill=n)) +
  theme_bw() +
  geom_tile(color='white') +
  scale_fill_viridis_c('',direction=-1) +
  labs(x='',y='') +
  coord_equal(ratio=0.8) +
  scale_y_continuous(breaks=seq(2002,2021,by=2))
# total of 13 sites
# in 2020-2021, 7 sites were visited over the year 
ggsave('Figures/Boysen_explore/Visits_moyr.png', height = 6.25, width=4.25, units='in')


### GIF! map with data frequency yearly ####
library(sf)
library(nhdplusTools)
start_point <- sf::st_as_sf(data.frame(x = -108.1847, y = 43.25167), 
                            coords = c("x", "y"), crs = 4326)
plot_nhdplus(start_point)


tmp <- Annualfreq |>
  filter(Latitude < 44)

Annualfreq.sf <- st_as_sf(tmp, coords=c('Longitude','Latitude'), crs=4269)


mapview::mapview(Annualfreq.sf)
ggplot(Annualfreq.sf) +
  geom_sf(aes(size=n))

# above is not working well
  
### Get NHD data! ####
library(rgdal)
library(raster)
## Flowlines
NHDflowline <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Spatial_Data/Boysen/NHD/NHDPLUS_H_1008_HU4_GDB.gdb', layer = 'NHDFlowline')
class(NHDflowline) # sf, df
crs(NHDflowline) # NAD83

## waterbody outlines
NHDwaterbody <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Spatial_Data/Boysen/NHD/NHDPLUS_H_1008_HU4_GDB.gdb', layer = 'NHDWaterbody')
class(NHDwaterbody)# sf, df
crs(NHDwaterbody) # NAD83

boysen <- st_crop(NHDwaterbody, xmin=-108.339946,ymin=43.142254,xmax=-108.039283,ymax=43.427151)
boysentribs <- st_crop(NHDflowline, xmin=-108.339946,ymin=43.142254,xmax=-108.039283,ymax=43.427151)

crs(Annualfreq.sf) == crs(boysen)


p<- ggplot() +
  geom_sf(boysentribs, mapping=aes(),color='cornflowerblue', alpha=0.5) +
  geom_sf(boysen, mapping=aes(),color='blue4', fill='cornflowerblue') +
  geom_sf(Annualfreq.sf |> filter(WaterbodyName != 'Muddy Creek Inlet'), 
          mapping=aes(size=n, group=`year(CollDate)`)) +
  theme_bw()

ggsave(p,'Figures/Bosyen_explore/alltimeMap.png', height=4.25, width=6.25, units='in')


# animating the map -- there is a bug in the package
library(gganimate)
library(transformr)
p
p.anim = p + transition_time(`year(CollDate)`)
animate(p.anim)
