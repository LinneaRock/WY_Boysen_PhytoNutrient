#----------------------------------#
# Script to call libraries and data 
#----------------------------------#

library(tidyverse)
library(lubridate)

# all WY chem data
ChemPhys <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                     fileEncoding="latin1") |>
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate)) |>
  rbind(read.csv('Data/RawData_WYDEQ/ChemPhysData_2022.csv', 
                 fileEncoding="latin1") |>
          mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
                 Year = year(CollDate)))

# filter for Boysen data
BoysenNutrient <- ChemPhys |>
  dplyr::select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  mutate(ChemUnits = sub('/l', '/L', ChemUnits)) |>
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  dplyr::select(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, BelowDet, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct() |>
  mutate(Year = year(CollDate)) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE))




BoysenChem <- ChemPhys |>
  dplyr::select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  mutate(ChemUnits = sub('/l', '/L', ChemUnits)) |>
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  dplyr::select(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, BelowDet, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(!ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct() |>
  mutate(Year = year(CollDate)) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE))



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
                 Year = year(CollDate)))


# filter for Boysen phyto data
BoysenPhyto <- BoysenChemPhys |>
  dplyr::select(StationID, WaterbodyName)  |>
  left_join(Phyto|>
              rename(chemName = WaterbodyName,
                     Genus.Species.Variety = `Genus/Species/Variety`)) |>
  drop_na(CollDate) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  unique()

