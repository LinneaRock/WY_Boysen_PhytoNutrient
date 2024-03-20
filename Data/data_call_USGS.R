
library(tidyverse)
library(dataRetrieval)

# outlet discharge downloaded from Bureau of Reclamation
 ## Observed Daily Average Lake/Reservoir Release - Total (cfs)
out_discharge <- read.csv('Data/boysen_outlet_discharge_BOR.csv', skip=7) |>
  mutate(CollDate = as.Date(Datetime..UTC.)) |>
  # average the twice daily measurements on some days
  group_by(CollDate) |>
  summarise(ChemValue = mean(Result)*28.316846592) |> # cfs to L/s
  ungroup() |>
  # rename things for consistency with DEQ data
  mutate(WaterbodyName = 'WIND RIVER BELOW BOYSEN RESERVOIR, WY',
         ShortName_Revised='Discharge',
         ChemUnits='l/s',
         Year=year(CollDate))
  


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
                           startDate = '1980-01-01')
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

# 00060 = discharge cfs
# 00010 = Water Temp degC
# 00095 = Spc uS/cm
# 00400 = pH
# 00600 = TN unfiltered mg/L
# 00605 = TON mg/L
# 00608 = Dissolved ammonia as N mg/L
# 00631 = Dissolved NO3 + NO2 as N mg/L
# 00665 = TP mg/L unfiltered
# 00671 = Dissolved orthophosphate as P mg/L
# 70507 = Total orthophosphate as P mg/L

site_ids <- c("06258000", "06253000", "06236100", '06259000')
parameterCd <- c('00095', '00400','00600', '00608', '00631', '00665','00671', '70507')

chemdata <- readWQPqw(paste0("USGS-", site_ids), parameterCd,
                     startDate = '1980-01-01')



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
                                     ShortName_Revised == 'Total Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)'~ 'Total Nitrogen (unfiltered)',
                                     ShortName_Revised == 'Dissolved Orthophosphate'~'Dissolved Orthophosphate',
                                     ShortName_Revised=='Discharge'~'Discharge',
                                     ShortName_Revised=='Dissolved Inorganic nitrogen (nitrate and nitrite)'~'Nitrate plus Nitrite as N',
                                     ShortName_Revised=='Dissolved Ammonia and ammonium'~'Total Ammonia as N')) |> # not actually total, but dissolved
  mutate(ChemUnits = case_when(ShortName_Revised=='Conductance'~'µS/cm',
                               ShortName_Revised=='Phosphorus as P (total)'~'mg/l',
                               ShortName_Revised=='pH'~'SU',
                               ShortName_Revised=='Total Nitrogen (unfiltered)'~'mg/l',
                               ShortName_Revised=='Dissolved Orthophosphate'~'mg/l',
                               ShortName_Revised=='Discharge'~'l/s',
                               ShortName_Revised=='Nitrate plus Nitrite as N'~'mg/l',
                               ShortName_Revised=='Total Ammonia as N'~'mg/l')) |>
  mutate(ChemValue=ifelse(ShortName_Revised=='Discharge', ChemValue*28.316846592, ChemValue)) |> # convert cfs to L/s 
bind_rows(out_discharge) |>
  mutate(Latitude=ifelse(WaterbodyName=='WIND RIVER BELOW BOYSEN RESERVOIR, WY',43.42496, Latitude),
         Longitude=ifelse(WaterbodyName=='WIND RIVER BELOW BOYSEN RESERVOIR, WY',-108.179, Longitude),
         DrainageArea=ifelse(WaterbodyName=='WIND RIVER BELOW BOYSEN RESERVOIR, WY', 7701, DrainageArea))

write.csv(BoysenTribs, 'Data/BoysenTribs.csv')

rm(list=c('BoysenTribs', 'chemdata', 'discharge','site_ids', 'parameterCd','windrvroutletparams','windrvrparams','fivemilecrkparams','muddycrkparams'))



