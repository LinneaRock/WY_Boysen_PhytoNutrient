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
  select(-X, -StationID)

reservoir <- ChemPhys |>
  filter(Year >= 2020) |>
  filter(grepl('Boysen', WaterbodyName)) |>
  filter(ShortName_Revised %in% c('Orthophosphate as P (total)', 'Total Ammonia as N', 'Phosphorus as P (total)', 'Total Nitrogen (unfiltered)', 'Nitrate plus Nitrite as N'))  |>
  mutate(ShortName_Revised = ifelse(ShortName_Revised=='Orthophosphate as P (total)', 'Orthophosphate', ShortName_Revised)) |>
  mutate(nutrient = ifelse(ShortName_Revised %in% c('Orthophosphate', 'Phosphorus as P (total)'), 'P', 'N')) |>
  mutate(type = ifelse(grepl('Transitional',WaterbodyName), 'Reservoir - transitional',
                       ifelse(grepl('Riverine',WaterbodyName), 'Reservoir - riverine',
                              ifelse(grepl('Lacustrine',WaterbodyName), 'Reservoir - lacustrine', 'Reservoir-bays and shores')))) |>
  select(-ChemSampID, -StationID) 

all_nutrient_data <- bind_rows(reservoir, tribs) |>
  select(-Comments, -SurfaceBottom, -MaximumDepth, -WBTypeID, -CommonName,-Section,-Quarter,-Range,-Town,-HUC12, -SampleDepth, -dec_coord_datum_cd, -Latitude, -Longitude, -DrainageArea, - Elevation, -WaterbodySize_Acres) |>
  mutate(BelowDet = ifelse(BelowDet ==1, 'Y', 'N')) |>
  filter(Year < 2023)

ggplot(all_nutrient_data |> filter(nutrient == 'P')) +
  geom_jitter(aes(CollDate, ChemValue, color = type),size=2) +
  geom_line(aes(CollDate, ChemValue,color = type),size=0.5) +
  scale_color_OkabeIto() +
  labs(x='', y='P concentration'~(mg~L^-1)) +
  facet_wrap(~ShortName_Revised, scales='free_y', ncol=1) +
  theme_classic() +
  theme(legend.title = element_blank())
ggsave('C:/Users/linne/OneDrive - University of Wyoming/Writing/proposal/ch3_P.png', height=4,width=6,dpi=1200)


ggplot(all_nutrient_data |> filter(nutrient == 'N')) +
  geom_jitter(aes(CollDate, ChemValue, color = type),size=2) +
  geom_line(aes(CollDate, ChemValue,color = type),size=0.5) +
  scale_color_OkabeIto() +
  labs(x='', y='N concentration'~(mg~L^-1)) +
  facet_wrap(~ShortName_Revised, scales='free_y', ncol=1) +
  theme_classic() +
  theme(legend.title = element_blank())
ggsave('C:/Users/linne/OneDrive - University of Wyoming/Writing/proposal/ch3_N.png', height=4,width=6,dpi=1200)




