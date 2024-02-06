#-----------------------------------------------#
# PCA attempts
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. Prep data for plotting ####
BoysenChemPhys <- ChemPhys |>
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
  mutate(year = year(CollDate)) |>
  filter(year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) 




BoysenPhyto <- BoysenChemPhys |>
  dplyr::select(StationID, WaterbodyName)  |>
  left_join(Phyto|>
              rename(chemName = WaterbodyName)) |>
  drop_na(CollDate) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  unique() |>
  group_by(StationID, WaterbodyName, CollDate, month, Year, `Genus/Species/Variety`) |>
  summarise(indsum = sum(`Individuals (Raw Cnt)`)) |>
  ungroup() |>
  group_by(StationID, WaterbodyName, month, Year) |>
  mutate(totaln=sum(indsum)) |>
  ungroup() |>
  mutate(perc =(indsum/totaln)*100) |>
  group_by(StationID, month, Year) |>
  mutate(check = sum(perc)) |>
  ungroup()

ggplot(BoysenPhyto) +
  geom_bar(aes(month, perc, fill=`Genus/Species/Variety`), stat='identity') +
  facet_wrap(WaterbodyName~Year)

n_phyto <- BoysenPhyto |>
  select(`Genus/Species/Variety`) |>
  distinct()

