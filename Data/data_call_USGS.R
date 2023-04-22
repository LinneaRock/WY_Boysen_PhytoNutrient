
library(tidyverse)
library(dataRetrieval)

# What data are available at each site? ####
muddycrkparams <- whatNWISdata(siteNumber = "06258000")
fivemilecrkparams <- whatNWISdata(siteNumber = "06253000")
windrvrparams <- whatNWISdata(siteNumber = "06236100")
windrvroutletparams <- whatNWISdata(siteNumber = '06259000')

# Load site data ####
sitedata <- readNWISsite(siteNumbers = c("06258000", "06253000", "06236100", '06259000'))
                         
# Load parameters we want ####
discharge <- readNWISdata(siteNumbers = c("06258000", "06253000", "06236100", '06259000'),
                           parameterCd = '00060',
                           startDate = '2002-01-01')
# chemdat <- readNWISqw(siteNumbers = c("06258000", "06253000", "06236100", '06259000'),
#                           parameterCd = c('00600', '00605', '00608', '00631', '00665','00671'),
#                           startDate = '2002-01-01') #|> readNWISqw will be removed. Need to use the code below to get water qual data from wqp:
#https://cran.r-project.org/web/packages/dataRetrieval/vignettes/qwdata_changes.html
  # mutate(parameter = case_when(param_cd == '00600' ~ 'TN_mgL', 
  #                         param_cd == '00605' ~ 'TON_mgL', 
  #                         param_cd == '00613' ~ 'Wind River Inlet', 
  #                         param_cd == '00631' ~ 'Wind River Outlet',
  #                         param_cd == '00665',
  #                         param_cd == '00671'))


site_ids <- c("06258000", "06253000", "06236100", '06259000')
parameterCd <- c('00010', '00095', '00400','00600', '00605', '00608', '00631', '00665','00671')

chemdata <- readWQPqw(paste0("USGS-", site_ids), parameterCd,
                     startDate = '2002-01-01')



# rename everything to match DEQ data
BoysenTribs <-discharge |>
              select(site_no, dateTime, X_00060_00003) |>
              mutate(ChemUnits='ft/s',
                     ShortName_Revised='Discharge')|>
              rename(ChemValue=X_00060_00003) |>
  bind_rows(chemdata |> 
              select(ActivityStartDate, MonitoringLocationIdentifier,CharacteristicName, ResultSampleFractionText, ResultMeasureValue, ResultMeasure.MeasureUnitCode) |>
              mutate(MonitoringLocationIdentifier = sub("^[^,]*-", "", MonitoringLocationIdentifier))|>
              unite(ShortName_Revised, ResultSampleFractionText, CharacteristicName, sep=' ') |>
              rename(dateTime=ActivityStartDate,
                     site_no=MonitoringLocationIdentifier,
                     ChemValue=ResultMeasureValue,
                     ChemUnits=ResultMeasure.MeasureUnitCode)) |>
  left_join(sitedata |>
               select(site_no, station_nm, dec_lat_va, dec_long_va, dec_coord_datum_cd, drain_area_va)) |>
  rename(StationID=site_no,
         WaterbodyName=station_nm,
         Latitude=dec_lat_va,
         Longitude=dec_long_va,
         DrainageArea=drain_area_va,
         CollDate=dateTime) |>
  group_by(StationID, CollDate, ChemUnits, ShortName_Revised, WaterbodyName, Latitude, Longitude, dec_coord_datum_cd, DrainageArea) |>
  summarise(ChemValue=mean(ChemValue)) |> # average any replicates
  ungroup() |>
  distinct() |>
  filter(ChemValue>=0) |>
  mutate(ShortName_Revised=case_when(ShortName_Revised =='Total Specific conductance' ~ 'Conductance',
                                     ShortName_Revised =='Total Phosphorus' ~ 'Phosphorus as P (total)',
                                     ShortName_Revised =='Total pH' ~ 'pH',
                                     ShortName_Revised == 'Total Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)', 'Total Nitrogen (unfiltered)')) |>
  filter(ShortName_Revised != 'Total Organic Nitrogen')

write.csv(BoysenTribs, 'Data/BoysenTribs.csv')

