#--------------------#
# Tributary loading 
#------------------- #

source('Data/CALL_DATA_LIB.R')


# 1. Fill in missing discharge data ####

ggplot(BoysenTribs |>
         filter(WaterbodyName %in% c('Muddy Creek', 'Wind River Inlet')) |>
         filter(ShortName_Revised == 'Discharge') |>
         mutate(fakedate = as.Date(paste0('1900-', month(CollDate), '-', day(CollDate))))) +
  geom_line(aes(fakedate, ChemValue, group=Year, color=Year)) +
  geom_point(BoysenTribs |>
               filter(WaterbodyName %in% c('Muddy Creek')) |>
               filter(ShortName_Revised == 'Discharge') |>
               filter(Year==1983) |>
               mutate(fakedate = as.Date(paste0('1900-', month(CollDate), '-', day(CollDate)))), 
             mapping=aes(fakedate, ChemValue)) + # 1983 weird outlier, do not include in average calculation below
  facet_wrap(~WaterbodyName, scales='free') +
  scale_color_viridis_c()

## 1a. Muddy Creek 2020 - missing Jan-Sep ####
mud <- BoysenTribs |>
  mutate(month=month(CollDate, label=TRUE, abbr=TRUE)) |>
  filter(WaterbodyName %in% c('Muddy Creek')) |>
  filter(ShortName_Revised == 'Discharge') |>
  filter(Year != 1983) |>
  filter(Year < 2024) |>
 # average data from 1980-82, 2020-23 (when we have discharge data, excluding 1983 super high discharge year)
  group_by(month) |>
  mutate(average_dis = mean(ChemValue)) |>
  select(-ChemValue, - CollDate, -Year) |>
  distinct() |>
  ungroup()


## 1b. Wind River Inlet 2020 - missing Jan-Sep ####
WRI <- BoysenTribs |>
  mutate(month=month(CollDate, label=TRUE, abbr=TRUE)) |>
  filter(WaterbodyName %in% c('Wind River Inlet')) |>
  filter(ShortName_Revised == 'Discharge') |>
  filter(Year < 2024) |>
  # average data from 1990-2013, 2020-23 
  group_by(month) |>
  mutate(average_dis = mean(ChemValue)) |>
  select(-ChemValue, - CollDate, -Year) |>
  distinct() |>
  ungroup()

## 1c. Fill in data ####

missingdatfill <- bind_rows(mud, WRI)


BoysenTribs_data <- BoysenTribs |>
  mutate(month=month(CollDate, label=TRUE, abbr=TRUE)) |>
  filter(between(Year, 2020,2023)) |> 
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN', 
                                       ShortName_Revised=='Total Ammonia as N'~'NH4', 
                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                       ShortName_Revised=='Dissolved Orthophosphate'~'PO4',
                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                       ShortName_Revised=='Discharge'~'Discharge')) |>
  drop_na(ShortName_Revised) |> # removing pH and SpC for now
  select(- ChemUnits) |> # all nutrients mg/L, Discharge L/s
  pivot_wider(names_from = ShortName_Revised, values_from = ChemValue) |>
  # monthly average discharge 
  group_by(WaterbodyName, Year, month) |>
  mutate(Discharge = mean(Discharge)) |>
  ungroup() |>
  filter(!if_all(c(TN, NH4, TP, PO4, NO3), is.na)) |>
  # fill in missing data with long-term averages
  left_join(missingdatfill) |>
  mutate(Discharge = ifelse(is.na(Discharge), average_dis, Discharge)) |>
  select(-average_dis, - ShortName_Revised, -ChemUnits, -CollDate) # all nutrients mg/L, Discharge L/s




# 2. Monthly mass loading, flux loading ####
# monthly time interval
q_interval <- 2.628e+6 # seconds/month

Trib_loading <- BoysenTribs_data |>
  # Mass loading kg/month: nutrint concentration (mg/L) * discharge (L/s) * interval(s/month) / convert mg to kg
  mutate(TN_load_kg = TN*Discharge*q_interval/1e6,
         TP_load_kg = TP*Discharge*q_interval/1e6,
         NH4_load_kg = NH4*Discharge*q_interval/1e6,
         NO3_load_kg = NO3*Discharge*q_interval/1e6,
         PO4_load_kg = PO4*Discharge*q_interval/1e6) |>
  select(-TN, -TP, -NH4, -NO3, -PO4, - Discharge) |>
  # Mass flux kg/ha/month: nutrint concentration (mg/L) * discharge (L/s) * interval(s/month) / convert mg to kg / DrainageArea (ha)
  mutate(TN_flux_kg_ha =  TN_load_kg/DrainageArea_ha,
         TP_flux_kg_ha =  TP_load_kg/DrainageArea_ha,
         NH4_flux_kg_ha = NH4_load_kg/DrainageArea_ha,
         NO3_flux_kg_ha = NO3_load_kg/DrainageArea_ha,
         PO4_flux_kg_ha = PO4_load_kg/DrainageArea_ha)

write.csv(Trib_loading, 'Data/trib_load_flux.csv')











