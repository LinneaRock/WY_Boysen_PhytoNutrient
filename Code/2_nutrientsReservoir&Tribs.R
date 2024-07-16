##########################################################
# Figure 2 nutrient loading and reservoir concentrations 
##########################################################

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')

# 2. tidy formatting of data for beautiful plotting ####
# load formatting 
Loading <- TribLoadFlux |> 
  mutate(IN_load_kg = NH4_load_kg + NO3_load_kg) |>
  select(-c(NH4_load_kg, NO3_load_kg, TotalTrib_NH4_kg, TotalTrib_NO3_kg, NH4_flux_kg_ha, NO3_flux_kg_ha)) |>
  mutate(fakedate = as.Date(paste0(Year,'-',month,'-01'), format='%Y-%b-%d')) |>
  distinct() |>
  pivot_longer(cols=c(TN_load_kg, TP_load_kg, IN_load_kg, PO4_load_kg), names_to = 'nutrient', values_to = 'load') |>
  mutate(nutrient=case_when(grepl('TN', nutrient)~'TN',
                            grepl('TP', nutrient)~'TP',
                            grepl('PO4', nutrient)~'Phosphate',
                            grepl('IN', nutrient)~'Inorganic N')) |>
  select(-c(7:12)) |>
  distinct() |>
  mutate(load=ifelse(grepl('Outlet', WaterbodyName), -1*load, load))

# flux formatting 
Fluxing <- TribLoadFlux |> 
  mutate(IN_flux_kg_ha = NH4_flux_kg_ha+NO3_flux_kg_ha) |>
  select(-c(NH4_load_kg, NO3_load_kg, TotalTrib_NH4_kg, TotalTrib_NO3_kg, NH4_flux_kg_ha, NO3_flux_kg_ha)) |>
  mutate(fakedate = as.Date(paste0(Year,'-',month,'-01'), format='%Y-%b-%d')) |>
  distinct() |>
  pivot_longer(cols=c(TN_flux_kg_ha, TP_flux_kg_ha, PO4_flux_kg_ha, IN_flux_kg_ha), names_to = 'nutrient', values_to = 'flux') |>
  mutate(nutrient=case_when(grepl('TN', nutrient)~'TN',
                            grepl('TP', nutrient)~'TP',
                            grepl('PO4', nutrient)~'Phosphate',
                            grepl('IN', nutrient)~'Inorganic N')) |>
  select(-c(7:12)) |>  
  distinct() |>
  mutate(flux=ifelse(grepl('Outlet', WaterbodyName), -1*flux, flux))


# total trib load formatting 
TotLoad <- TribLoadFlux |> 
  mutate(TotalTrib_IN_kg = TotalTrib_NH4_kg + TotalTrib_NO3_kg) |> 
  select(-c(NH4_load_kg, NO3_load_kg, TotalTrib_NH4_kg, TotalTrib_NO3_kg, NH4_flux_kg_ha, NO3_flux_kg_ha)) |>
  mutate(fakedate = as.Date(paste0(Year,'-',month,'-01'), format='%Y-%b-%d')) |>
  distinct() |>
  pivot_longer(cols=c(TotalTrib_IN_kg, TotalTrib_TN_kg, TotalTrib_TP_kg, TotalTrib_PO4_kg), names_to = 'nutrient', values_to = 'TotalLoad') |>
  mutate(nutrient=case_when(grepl('TN', nutrient)~'TN',
                            grepl('TP', nutrient)~'TP',
                            grepl('PO4', nutrient)~'Phosphate',
                            grepl('IN', nutrient)~'Inorganic N')) |>
  select(-c(6:12))
  

# combine data
all_load_data <- left_join(Loading, Fluxing) |>
  left_join(TotLoad) |>
  distinct() |>
  mutate(dummy = ifelse(is.na(load), 'a','b')) |>
  arrange(dummy,fakedate,WaterbodyName) |>
  mutate(nutrient=factor(nutrient, levels=c('TN','TP', 'Inorganic N', 'Phosphate')))
  

# 3. Plot tribs! ####
trib_colors <- c('#117733','#DDCC77','#882255','#332288')

ggplot(all_load_data) +
  geom_point(aes(fakedate, load, color=WaterbodyName)) +
  geom_path(aes(fakedate, load, color=WaterbodyName)) +
  scale_color_manual('',values =trib_colors) +
  geom_point(aes(fakedate, TotalLoad),shape=21) +
  facet_wrap(~nutrient, scales='free_y') +
  labs(x='',y='Monthly nutrient load (kg)') +
  theme_minimal() +
  theme(legend.position = 'bottom')
  
  
# 4. What percent of nutrients into reservoir are retained? ####
retention <- all_load_data |>
  select(-c(1:6, 10, 12))
  

