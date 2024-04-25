#----------------------------------#
# Script to call libraries and data 
#----------------------------------#

library(tidyverse)
library(lubridate)
library(patchwork)
library(readxl)

# call in schmidt's stability 
SS <- read.csv('Data/Schmidts_Stability.csv') |>
  select(-X) |>
  mutate(CollDate = as.Date(CollDate)) |>
  mutate(ChemUnits='J_m2') |>
  select(-Year, -Month)


# all WY chem data
ChemPhys <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                     fileEncoding="latin1") |>
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate)) |>
  rbind(read.csv('Data/RawData_WYDEQ/ChemPhysData_2022.csv', 
                 fileEncoding="latin1") |>
          mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
                 Year = year(CollDate))) |>
  bind_rows(read_xlsx('Data/RawData_WYDEQ/2023_data/Boysen_2023_Chemistry_Phytoplankton.xlsx', sheet=1) |>
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate),
         HUC12=as.numeric(HUC12)) ) |>
  # for all values below detection, make them half the reporting limit. In all cases of below detection, the reported value is the reporting limit. It is not always consistent, sometimes due to sample volume limits, matrix interference, laboratory instrument limitations, etc., the limits for individual analyte results are more or less than standard reporting limits. <- from personal communication with chloe Schaub and Eric Hargett at DEQ April 2024
  mutate(ChemValue = ifelse(BelowDet==1, ChemValue/2,ChemValue)) |>
  mutate(Year = year(CollDate)) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  mutate(julianday=yday(CollDate))

# filter for Boysen data
BoysenNutrient <- ChemPhys |>
  dplyr::select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  mutate(ChemUnits = sub('/l', '/L', ChemUnits)) |>
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  dplyr::select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, ShortName_Revised, BelowDet, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N", "Orthophosphate as P (dissolved)")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct() 


BoysenNutrient$WaterbodyName = trimws(BoysenNutrient$WaterbodyName) 

SS <- SS |>
  left_join(BoysenNutrient |> select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday)) |>
  distinct()



BoysenChem <- ChemPhys |>
  dplyr::select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  mutate(ChemUnits = sub('/l', '/L', ChemUnits)) |>
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  dplyr::select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, ShortName_Revised, BelowDet, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(!ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N", "Orthophosphate as P (dissolved)")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct()  

BoysenChem$WaterbodyName = trimws(BoysenChem$WaterbodyName)


BoysenChem <- BoysenChem |> 
  bind_rows(SS |> select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, SS, ChemUnits) |> rename(ChemValue=SS) |> mutate(ShortName_Revised='Stability'))




# tributary and outlet data
BoysenTribs <- read.csv('Data/BoysenTribs.csv') |>
  mutate(CollDate = as.Date(CollDate, format='%Y-%m-%d'),
         Year = year(CollDate)) |>
  select(-X, -dec_coord_datum_cd, -StationID) |>
  mutate(WaterbodyName = case_when(WaterbodyName=='WIND RIVER BELOW BOYSEN RESERVOIR, WY'~'Wind River Outlet',
                                   WaterbodyName=='WIND RIVER AB BOYSEN RESERVOIR, NR SHOSHONI, WY'~'Wind River Inlet',
                                   WaterbodyName=='FIVEMILE CREEK NEAR SHOSHONI, WY'~'Fivemile Creek',
                                   WaterbodyName=='MUDDY CREEK NEAR SHOSHONI, WY'~'Muddy Creek'))


library(readr)
# all WY phyto data
Phyto <- read_csv("Data/RawData_WYDEQ/Phytoplankton_2013-2021.csv", 
                                    skip = 5) |>
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate)) |>
  rbind(read_csv("Data/RawData_WYDEQ/Phytoplankton_2022.csv", 
                 skip = 5) |>
          mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
                 Year = year(CollDate))) |>
  filter(Division != 'Copepoda',
                   !`Genus/Species/Variety` %in% c('Daphnia','Cladocera', '	
          Copepoda: nauplius')) |>
  bind_rows(read_xlsx('Data/RawData_WYDEQ/2023_data/Boysen_2023_Chemistry_Phytoplankton.xlsx', sheet=2) |>
          mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
                 Year = year(CollDate)))


# filter for Boysen phyto data
BoysenPhyto <- BoysenNutrient |>
  dplyr::select(StationID, WaterbodyName)  |>
  left_join(Phyto|>
              rename(chemName = WaterbodyName,
                     Genus.Species.Variety = `Genus/Species/Variety`)) |>
  drop_na(CollDate) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  unique()




# get profile data and do annoying things to get names right
BoysenProfile <- read.csv('Data/profiles.csv') |>
  mutate(CollDate=as.Date(CollDate, format='%m/%d/%Y'))
#there's a space at the beginning of all the waterbody names :( 
annoyingworkaroundfornames <- sub('.', '', BoysenProfile$WaterbodyName)
BoysenProfile <- BoysenProfile |>
  select(-WaterbodyName) |>
  mutate(WaterbodyName=annoyingworkaroundfornames)
  
rm(annoyingworkaroundfornames)







