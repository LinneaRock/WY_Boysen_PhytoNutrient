# proposal figures 

source('Data/CALL_DATA_LIB.R')
library(colorblindr)

tribs <- BoysenTribs |> 
  filter(Year >= 2020) |>
  mutate(type = ifelse(WaterbodyName == 'WIND RIVER BELOW BOYSEN RESERVOIR, WY', 'Outlet', 
                       ifelse(WaterbodyName == 'WIND RIVER AB BOYSEN RESERVOIR, NR SHOSHONI, WY', 'Main inlet', 'Other inlets'))) |>
  filter(ShortName_Revised %in% c('Dissolved Orthophosphate','Phosphorus as P (total)', 'Total Nitrogen (unfiltered)')) |>
  mutate(ShortName_Revised = ifelse(ShortName_Revised=='Dissolved Orthophosphate', 'Orthophosphate', ShortName_Revised)) |>
  mutate(nutrient = ifelse(ShortName_Revised %in% c('Orthophosphate', 'Phosphorus as P (total)'), 'P', 'N')) |>
  select(-X, -StationID) |>
  mutate(ecotype = 'Rivers')


reservoir <- ChemPhys |>
  filter(Year >= 2020) |>
  filter(grepl('Boysen', WaterbodyName)) |>
  filter(ShortName_Revised %in% c('Orthophosphate as P (total)', 'Total Ammonia as N', 'Phosphorus as P (total)', 'Total Nitrogen (unfiltered)', 'Nitrate plus Nitrite as N'))  |>
  mutate(ShortName_Revised = ifelse(ShortName_Revised=='Orthophosphate as P (total)', 'Orthophosphate', ShortName_Revised)) |>
  mutate(nutrient = ifelse(ShortName_Revised %in% c('Orthophosphate', 'Phosphorus as P (total)'), 'P', 'N')) |>
  mutate(type = ifelse(grepl('Transitional',WaterbodyName), 'Reservoir - transitional',
                       ifelse(grepl('Riverine',WaterbodyName), 'Reservoir - riverine',
                              ifelse(grepl('Lacustrine',WaterbodyName), 'Reservoir - lacustrine', 'Reservoir-bays and shores')))) |>
  select(-ChemSampID, -StationID) |>
  mutate(ecotype = 'In-reservoir') |>
  mutate(detectionlimit = ifelse(ShortName_Revised == 'Total Ammonia as N', 0.05,
                                 ifelse(ShortName_Revised ==  'Total Nitrogen (unfiltered)', 0.1,
                                        ifelse(ShortName_Revised == 'Nitrate plus Nitrite as N', 0.05,
                                               ifelse(ShortName_Revised == 'Orthophosphate', 0.01,
                                                      ifelse(ShortName_Revised == 'Phosphorus as P (total)', 0.01, NA))))))

#na_dates <- data.frame(CollDate = seq.Date(as.Date('2020-01-01'), as.Date('2022-12-01'), by='month'))

boysensites <- unique(reservoir$WaterbodyName)
boysentribs <- unique(BoysenTribs$WaterbodyName)


all_nutrient_data <- bind_rows(reservoir, tribs) |>
  select(-Comments, -SurfaceBottom, -MaximumDepth, -WBTypeID, -CommonName,-Section,-Quarter,-Range,-Town,-HUC12, -SampleDepth, -dec_coord_datum_cd, -Latitude, -Longitude, -DrainageArea, - Elevation, -WaterbodySize_Acres) |>
  mutate(BelowDet = ifelse(BelowDet ==1, 'Y', 'N')) |>
  filter(Year < 2023) |>
  filter(month(CollDate) %in% c(5:10))
 # mutate(monyear = as.Date(paste0(year(CollDate), '-', month(CollDate), '-', '1')))
  # pivot_wider(id_cols = c(ShortName_Revised, CollDate, ChemUnits, nutrient), names_from = 'WaterbodyName', values_from = 'ChemValue') |>
  # pivot_longer(5:15, names_to = 'WaterbodyName', values_to = 'ChemValue') |>
  # mutate(ecotype = ifelse(WaterbodyName %in% boysentribs, 'Rivers', 'In-reservoir')) |>
  # mutate(type = ifelse(WaterbodyName %in% boysensites & grepl('Transitional',WaterbodyName), 'Reservoir - transitional',
  #                      ifelse(WaterbodyName %in% boysensites &grepl('Riverine',WaterbodyName), 'Reservoir - riverine',
  #                             ifelse(WaterbodyName %in% boysensites & grepl('Lacustrine',WaterbodyName), 'Reservoir - lacustrine', 'Reservoir-bays and shores')))) |>
  # mutate(type = ifelse(WaterbodyName == 'WIND RIVER BELOW BOYSEN RESERVOIR, WY', 'Outlet', 
  #                      ifelse(WaterbodyName == 'WIND RIVER AB BOYSEN RESERVOIR, NR SHOSHONI, WY', 'Main inlet', 
  #                             ifelse(WaterbodyName %in% c('MUDDY CREEK NEAR SHOSHONI, WY', 'FIVEMILE CREEK NEAR SHOSHONI, WY'), 'Other inlets', type)))) |>
  # mutate(fakeDate = as.Date(paste(year(CollDate), month(CollDate), '01', sep = '-'))) |>
  # select(-CollDate) |>
  # distinct()

ggplot(all_nutrient_data |> filter(nutrient == 'P')) +
  geom_jitter(aes(CollDate, ChemValue, color = type),size=2) +
  #geom_line(aes(CollDate, ChemValue, color = type),size=0.5) +
  scale_color_OkabeIto() +
  labs(x='', y='P concentration'~(mg~L^-1)) +
  facet_wrap(~ShortName_Revised*ecotype, scales='free', ncol=2) +
  geom_line(aes(CollDate, detectionlimit)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_x_date(date_labels = '%Y') +
  annotate("rect", xmin = as.Date('2020-11-01'), xmax = as.Date('2021-05-01'), ymin = 0, ymax = Inf, alpha=0.2, color = "grey") +
  annotate("rect", xmin = as.Date('2021-11-01'), xmax = as.Date('2022-05-01'), ymin = 0, ymax = Inf, alpha=0.2, color = "grey") 
 ggsave('C:/Users/linne/OneDrive - University of Wyoming/Writing/proposal/ch3_P.png', height=6,width=8,dpi=1200)



ggplot(all_nutrient_data |> filter(nutrient == 'N')) +
  geom_jitter(aes(CollDate, ChemValue, color = type),size=2) +
  #geom_line(aes(CollDate, ChemValue,color = type),size=0.5) +
  scale_color_OkabeIto() +
  labs(x='', y='N concentration'~(mg~L^-1)) +
  facet_wrap(~ShortName_Revised*ecotype, scales='free', ncol=2) +
  geom_line(aes(CollDate, detectionlimit)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_x_date(date_labels = '%Y') +
  annotate("rect", xmin = as.Date('2020-11-01'), xmax = as.Date('2021-05-01'), ymin=0, ymax = Inf, alpha=0.2, color = "grey") +
  annotate("rect", xmin = as.Date('2021-11-01'), xmax = as.Date('2022-05-01'), ymin = 0, ymax = Inf, alpha=0.2, color = "grey") 
ggsave('C:/Users/linne/OneDrive - University of Wyoming/Writing/proposal/ch3_N.png', height=6,width=8,dpi=1200)






