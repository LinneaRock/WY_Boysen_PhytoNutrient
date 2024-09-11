#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Tributary nutrient loading and reservoir concentrations 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')

# TRIBUTARY LOADING ####
## 2. tidy formatting of data for beautiful plotting ####
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
  

## 3. Plot tribs loading! ####
options(scipen = 999)

trib_colors <- c('#117733','#DDCC77','#882255','#332288')

ggplot(all_load_data) +
  geom_point(aes(fakedate, load, color=WaterbodyName)) +
  geom_path(aes(fakedate, load, color=WaterbodyName)) +
  scale_color_manual('',values =trib_colors) +
  geom_point(aes(fakedate, TotalLoad),shape=21) +
  facet_wrap(~nutrient, scales='free_y') +
  labs(x='',y='Nutrient mass (kg)') +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave('Figures/nutrientdynamics_tribs.png', height=4.5, width=6.5, units='in',dpi=1200)

### 3a. discharge fig ####
discharge <- BoysenTribs |>
  filter(between(year(CollDate), 2020, 2023)) |>
  filter(ShortName_Revised=='Discharge') |>
  mutate(ChemValue=ifelse(WaterbodyName=='Wind River Outlet', ChemValue*-1, ChemValue))

estimated_dis <- estimated_discharge |> 
  filter(estimated=='Y') |>
  mutate(CollDate = as.Date(paste(month, '15', Year, sep='-'), format='%b-%d-%Y'))

ggplot() + 
  geom_line(discharge, mapping=aes(CollDate, ChemValue, color=WaterbodyName)) +
  geom_line(estimated_dis, mapping=aes(CollDate, Discharge, color=WaterbodyName), linetype=2) +
  scale_color_manual('',values =trib_colors) +
  labs(x='',y='Discharge '~(L~s^-1)) +
  theme_minimal() +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.1, 0.85)) +
  annotate('text', x=as.Date('2020-01-01'), y= 200000, label='- - Estimated discharge', hjust=0.15, size=3.3)
ggsave('Figures/discharge.png', height=4.5, width=6.5, units='in',dpi=1200)


### 3a. nutrient concentration fig ####
trib_nuts <- BoysenTribs_data |>
  mutate(`Inorganic N`=NO3 + NH4) |>
  rename(Phosphate=PO4) |>
  mutate(CollDate = as.Date(paste(month, '15', Year, sep='-'), format='%b-%d-%Y')) |>
  pivot_longer(cols=c('TN', `Inorganic N`, 'TP', 'Phosphate')) |>
  mutate(name = factor(name, levels=c('TN','TP','Inorganic N', 'Phosphate')))

ggplot(trib_nuts) + 
  geom_point(aes(CollDate, value, color=WaterbodyName)) +
  geom_line(aes(CollDate, value, color=WaterbodyName)) +
  scale_color_manual('',values =trib_colors) +
  facet_wrap(~name, scales = 'free') +
  labs(x='',y='Concentration '~(mg~L^-1)) +
  theme_minimal() 
ggsave('Figures/tribs_nutrientConcentrations.png', height=4.5, width=6.5, units='in',dpi=1200)

## 4. Plot area-normalized flux differences ####
ggplot(all_load_data) +
  geom_boxplot(aes(WaterbodyName, flux, color=WaterbodyName)) +
  scale_color_manual('',values =trib_colors) +
  facet_wrap(~nutrient, scales='free_y')
  


## 5. What percent of nutrients into reservoir are retained? ####
tmp <- all_load_data |>
  #select(-c(1:6, 10, 12)) |>
  mutate(load=ifelse(!is.na(TotalLoad), NA, load)) |>
  distinct() 

tot.tmp <- tmp |>
  select(TotalLoad, fakedate, nutrient) |>
  distinct() |>
  drop_na()

out.tmp <- tmp |>
  select(load, fakedate, nutrient) |>
  filter(load<0) |>
  distinct() |>
  drop_na()


retention <- left_join(tot.tmp, out.tmp) |> # we will only report on data we can actual calculate something
  mutate(subtraction = TotalLoad + load) |> # how much entering reservoir does not exit
  mutate(perc_ret = ifelse(subtraction > 0, abs(load)/TotalLoad * 100, load/TotalLoad *100)) |>
  group_by(nutrient) |>
  mutate(median_perc_ret = median(perc_ret, na.rm=TRUE)) |>
  ungroup() |>
  mutate(m = month(fakedate)) |>
  mutate(month=month(fakedate, label=TRUE))


summary(lm(perc_ret~m, retention |> filter(nutrient=='TN')))
#slope=17.52
#p-value: 0.0007671
#Adjusted R-squared:  0.6267 

summary(lm(perc_ret~m, retention |> filter(nutrient=='TP')))
#slope=-8.924
#p-value: 0.4045
#Adjusted R-squared:  -0.00957 

summary(lm(perc_ret~m, retention |> filter(nutrient=='Inorganic N')))
#slope=0.1973 
#p-value: 0.9614
#Adjusted R-squared:  -0.4978 

summary(lm(perc_ret~m, retention |> filter(nutrient=='Phosphate')))
#slope=-10.057
#p-value: 0.5918
#Adjusted R-squared:  -0.04298 


ggplot(retention) +
  geom_boxplot(aes(month, perc_ret, group=month)) +
  geom_smooth(aes(m, perc_ret), method='lm', se=FALSE) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  facet_wrap(~nutrient) +
  labs(x='',y='Nutrient retention (%)')
ggsave('Figures/retentionTimeFig.png', height=4.5, width=6.5, units='in',dpi=1200)




# IN-RESERVOIR NUTRIENT CONCENTRATIONS ####
## 6. Tidy data for beautiful plotting ####
reservoir_nutrients <- BoysenNutrient |>
  filter(ShortName_Revised %in% c('TN', 'TP', 'NO3', 'NH4', 'PO4')) |>
  pivot_wider(names_from = ShortName_Revised, values_from='ChemValue') |>
  mutate(`Inorganic N` = NH4 + NO3) |>
  select(-NH4, -NO3) |>
  rename(Phosphate=PO4) |>
  pivot_longer(cols=c(Phosphate, TP, TN, `Inorganic N`), names_to = 'nutrient', values_to = 'concentration') |>
  group_by(CollDate, nutrient) |>
  mutate(meanConc = mean(concentration)) |>
  ungroup() |>
  mutate(nutrient=factor(nutrient, levels=c('TN','TP', 'Inorganic N', 'Phosphate'))) |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay')))


## 7. Plot in-reservoir nutrients ####
ggplot(reservoir_nutrients) +
  geom_point(aes(CollDate, concentration, color=WaterbodyName), alpha=0.5) +
  #geom_path(aes(CollDate, concentration, color=WaterbodyName), alpha=0.5) +
  scale_color_viridis_d('', option='turbo') +
  geom_point(aes(CollDate, meanConc),shape=22, size=2) +
  geom_line(aes(CollDate, meanConc, group=Year)) +
  facet_wrap(~nutrient, scales='free_y') +
  labs(x='',y='Nutrient concentration'~(mg~L^-1)) +
  theme_minimal() 
ggsave('Figures/nutrientdynamics_boysen.png', height=4.5, width=6.5, units='in',dpi=1200)

# # 8. Put plots together ####
# a/b +
#   plot_annotation(tag_levels = 'a', tag_suffix = ')')
# ggsave('Figures/nutrientdynamics.png', height=10.5, width=8.5, units='in',dpi=1200)
  
