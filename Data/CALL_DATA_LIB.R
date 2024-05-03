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



# get profile data and do annoying things to get names right
BoysenProfile <- read.csv('Data/profiles.csv') |>
  mutate(CollDate=as.Date(CollDate, format='%m/%d/%Y'))
#there's a space at the beginning of all the waterbody names :( 
annoyingworkaroundfornames <- sub('.', '', BoysenProfile$WaterbodyName)
BoysenProfile <- BoysenProfile |>
  select(-WaterbodyName) |>
  mutate(WaterbodyName=annoyingworkaroundfornames) |>
  group_by(WaterbodyName, CollDate) |>
  mutate(maxdepth = max(depth_m)) |>
  ungroup()

rm(annoyingworkaroundfornames)


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
  select(-BelowDet) |>
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
  dplyr::select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, ShortName_Revised, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N", "Orthophosphate as P (dissolved)")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct() |>
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN',
                                       ShortName_Revised=='Total Ammonia as N'~'NH4',
                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                       ShortName_Revised=='Orthophosphate as P (total)'~'PO4',
                                       ShortName_Revised=='Orthophosphate as P (dissolved)'~'PO4',
                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                       ShortName_Revised=='Chlorophyll a (phytoplankton)'~'CHLA')) 


stoich <- BoysenNutrient |>
  select(-ChemUnits,-SampleDepth) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) |>
  mutate(TN.TP = (TN/TP)*2.11306,
         IN.PO4 = ((NH4+NO3)/PO4)*2.11306) |>
  select(-c(PO4,NH4,TP,TN,NO3,CHLA)) |>
  mutate(ChemUnits='molar ratio') |>
  pivot_longer(c(TN.TP, IN.PO4), names_to = 'ShortName_Revised', values_to = 'ChemValue')

BoysenNutrient <- BoysenNutrient |>
  bind_rows(stoich)

rm(stoich)
  


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
  dplyr::select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, ShortName_Revised, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(!ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N", "Orthophosphate as P (dissolved)")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct() |>
  mutate(ShortName_Revised = case_when(
    ShortName_Revised == 'Secchi Depth' ~'Secchi',
    ShortName_Revised=='pH'~'pH',
    ShortName_Revised=='DO, mg/L'~'DO',
    ShortName_Revised=='Conductance'~'SpC',
    ShortName_Revised=='Stability'~'Stability',
    ShortName_Revised=='Temp'~'Temp')) |>
  # getting rid of other variables for now
  drop_na(ShortName_Revised)

BoysenChem$WaterbodyName = trimws(BoysenChem$WaterbodyName)


BoysenChem <- BoysenChem |> 
  bind_rows(SS |> select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, SS, ChemUnits) |> rename(ChemValue=SS) |> mutate(ShortName_Revised='Stability')) |>
  bind_rows(BoysenProfile |> left_join(BoysenNutrient |> select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday)) |> select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday, maxdepth) |> distinct() |> pivot_longer(maxdepth, names_to = 'ShortName_Revised', values_to = 'ChemValue') |> mutate(ChemUnits='m'))



# remove funky data from DO and temp

BoysenChem <- BoysenChem |>
  #Lacustrine Pelagic: Dam 2021-09-15 DO doesn't make sense
  mutate(ChemValue=ifelse(ShortName_Revised=='DO' & ChemValue<0.5, NA, ChemValue), # looking at the profiles, the DO values are near zero through water column this day. seems like something is wrong 
    #Tough Creek Campground 2023-08-16 - this was definitely a misreport or faulty sensor. 
         ChemValue=ifelse(ShortName_Revised=='Temp' & ChemValue<6, NA, ChemValue))

# replace funky data manually:
## from profile for TCC temperature
replace_value_temp <- as.numeric(BoysenProfile |>
  filter(WaterbodyName =='Tough Creek Campground' &
           CollDate==as.Date('2023-08-16') &
           between(depth_m, 0, 1)) |>
  summarise(mean(temp_C)))

## using average value from other data collected at surface in September in LPD
replace_value_do <- as.numeric(BoysenChem |>
                                   filter(WaterbodyName =='Lacustrine Pelagic: Dam' &
                                            month=='Sep',
                                          ShortName_Revised=='DO') |>
                                   summarise(mean(ChemValue, na.rm = TRUE)))

BoysenChem <- BoysenChem |>
  mutate(ChemValue= ifelse(ShortName_Revised=='Temp' & is.na(ChemValue), replace_value_temp, ChemValue),
         ChemValue= ifelse(ShortName_Revised=='DO' & is.na(ChemValue), replace_value_do, ChemValue))

rm(replace_value_do)
rm(replace_value_temp)

# tributary and outlet data
BoysenTribs <- read.csv('Data/BoysenTribs.csv') |>
  mutate(CollDate = as.Date(CollDate, format='%Y-%m-%d'),
         Year = year(CollDate)) |>
  select(-X, -dec_coord_datum_cd, -StationID) |>
  mutate(WaterbodyName = case_when(WaterbodyName=='WIND RIVER BELOW BOYSEN RESERVOIR, WY'~'Wind River Outlet',
                                   WaterbodyName=='WIND RIVER AB BOYSEN RESERVOIR, NR SHOSHONI, WY'~'Wind River Inlet',
                                   WaterbodyName=='FIVEMILE CREEK NEAR SHOSHONI, WY'~'Fivemile Creek',
                                   WaterbodyName=='MUDDY CREEK NEAR SHOSHONI, WY'~'Muddy Creek')) |>
  mutate(DrainageArea_ha = DrainageArea * 258.999) |>
  select(-DrainageArea)

TribLoadFlux <- read.csv('Data/trib_load_flux.csv') |>
  select(-X)


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




# read in cyano toxin data and filter for presence/absence 
cyanotoxin <- read_xlsx('Data/cyano_data.xlsx') |>
  mutate(CollDate=as.Date(CollDate),
         Anatoxin_a_ugL=as.numeric(Anatoxin_a_ugL),
         Microcysitn_ugL=as.numeric(Microcysitn_ugL)) |>
  mutate(toxinpresent = ifelse(is.na(Anatoxin_a_ugL) & is.na(Microcysitn_ugL), 'N', 'Y')) |>
  select(-Anatoxin_a_ugL, -Microcysitn_ugL) |>
  filter(toxinpresent=='Y' |
           !is.na(Advisory)) |>
  distinct() 



# combine data for easier to use all boysen dataframe
phyto_class <- read.csv('Data/phyto_class.csv') 

BoysenPhyto_A <- BoysenPhyto |>
  left_join(phyto_class)  |>
  filter(RepNum == 0) |> #only keep first counts  to avoid the insane confusion I had when I ignored this column
  distinct()|>
  group_by(WaterbodyName, CollDate, month, Year, Genus.Species.Variety) |>
  mutate(indsum = sum(`Individuals (Raw Cnt)`),
         biovolumeSum_cellsL=sum(`Density (cells/L)`)) |>
  ungroup() |>
  distinct() |>
  group_by(WaterbodyName, month, Year) |>
  mutate(totaln=sum(indsum),
         totalbiovolume_cellsL=sum(biovolumeSum_cellsL)) |>
  ungroup() |>
  distinct() |>
  mutate(percCount =(indsum/totaln)*100,
         percBiovolume=(biovolumeSum_cellsL/totalbiovolume_cellsL)*100) |>
  group_by(WaterbodyName, month, Year) |>
  mutate(checkn = sum(percCount),
         checkbv=sum(percBiovolume)) |>
  ungroup() |>
  select(WaterbodyName, CollDate, month, Year, Genus.Species.Variety, indsum, totaln, percCount, checkn,biovolumeSum_cellsL, totalbiovolume_cellsL,percBiovolume,checkbv) |>
  # percent biovolume = percent count -- duh but needed to verify
  
  # Shannon-Weiner Diversity Index; <1.5 low diversity, >2.5 high diversity'
  group_by(WaterbodyName, CollDate) |>
  mutate(H=sum(abs(log(indsum/totaln)*(indsum/totaln)))) |>
  ungroup()


H <- BoysenPhyto_A |>
  select(WaterbodyName,CollDate,H) |>
  left_join(BoysenNutrient |> select(StationID, WaterbodyName, Latitude, Longitude, CollDate, Year, month, julianday)) |>
  distinct() |>
  pivot_longer(H, names_to = 'ShortName_Revised',values_to = 'ChemValue') 


BoysenChem <- BoysenChem |>
  bind_rows(H)


BoysenPhyto_cat <- BoysenPhyto |>
  filter(RepNum == 0) |>
  left_join(phyto_class)  |>
  left_join(BoysenPhyto_A |> select(WaterbodyName, Year, H, month, totaln, totalbiovolume_cellsL)) |>
  distinct() |>
  group_by(WaterbodyName, CollDate, month, Year, cat, totaln, totalbiovolume_cellsL, H) |>
  summarise(catsum = sum(`Individuals (Raw Cnt)`),
            CATbiovolumeSum_cellsL=sum(`Density (cells/L)`)) |>
  ungroup() |>
  mutate(Catperc =(catsum/totaln)*100,
         CATpercBiovolume=(CATbiovolumeSum_cellsL/totalbiovolume_cellsL)*100) |>
  group_by(WaterbodyName, month, Year,H) |>
  mutate(checkCat = sum(Catperc),
         checkbv=sum(CATpercBiovolume)) |>
  # percent biovolume = percent count -- duh but needed to verify
  ungroup() |>
  select(WaterbodyName, month, Year, H, cat, Catperc) |>
  distinct() |>
  pivot_wider(names_from = cat, values_from = Catperc) 
BoysenPhyto_cat <- replace(BoysenPhyto_cat, is.na(BoysenPhyto_cat), 0)


rm(H)
