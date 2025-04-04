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
  mutate(load=ifelse(grepl('Outlet', WaterbodyName), -1*load, load)) |>
  distinct() |>
  mutate(WaterbodyName = factor(WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet')))

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
  mutate(flux=ifelse(grepl('Outlet', WaterbodyName), -1*flux, flux)) |>
  mutate(WaterbodyName = factor(WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet')))


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
  select(-c(6:12)) |>
  mutate(WaterbodyName = factor(WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet')))
  

# combine data
all_load_data <- left_join(Loading, Fluxing) |>
  left_join(TotLoad) |>
  distinct() |>
  mutate(dummy = ifelse(is.na(load), 'a','b')) |>
  arrange(dummy,fakedate,WaterbodyName) |>
  mutate(nutrient=factor(nutrient, levels=c('TN','TP', 'Inorganic N', 'Phosphate'))) |>
  mutate(WaterbodyName = factor(WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet')))
  

## 3. Plot tribs loading! ####
options(scipen = 999)

trib_colors <- c('#332288','#882255','#E6AA68','#DDCC77')

ggplot(all_load_data) +
  geom_point(aes(fakedate, load, color=WaterbodyName)) +
  geom_path(aes(fakedate, load, color=WaterbodyName)) +
  scale_color_manual('',values =trib_colors) +
  geom_point(aes(fakedate, TotalLoad),shape=21, size=2.5) +
  facet_wrap(~nutrient, scales='free_y') +
  labs(x='',y='Nutrient mass (kg)') +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave('Figures/nutrientdynamics_tribs.png', height=4.5, width=6.5, units='in',dpi=1200)

### 3a. discharge fig ####
discharge <- BoysenTribs |>
  filter(between(year(CollDate), 2020, 2023)) |>
  filter(ShortName_Revised=='Discharge') |>
  mutate(ChemValue=ifelse(WaterbodyName=='Wind River Outlet', ChemValue*-1, ChemValue)) |>
  mutate(WaterbodyName = factor(WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet')))

estimated_dis <- estimated_discharge |> 
  filter(estimated=='Y') |>
  mutate(CollDate = as.Date(paste(month, '15', Year, sep='-'), format='%b-%d-%Y'))

ggplot() + 
  geom_line(discharge, mapping=aes(CollDate, ChemValue, color=WaterbodyName)) +
  geom_line(estimated_dis, mapping=aes(CollDate, Discharge, color=WaterbodyName), linetype=2) +
  scale_color_manual('',values =trib_colors) +
  labs(x='',y='Discharge '~(L~s^-1)) +
  theme_minimal() +
  theme(legend.position = 'none') +
  annotate('text', x=as.Date('2020-01-01'), y= 250000, label='- - Estimated discharge', hjust=0.15, size=3.3)
ggsave('Figures/discharge.png', height=4.5, width=6.5, units='in',dpi=1200)


### 3a. nutrient concentration fig ####
trib_nuts <- TribLoadFlux |>
  mutate(`Inorganic N`=NO3 + NH4) |>
  rename(Phosphate=PO4) |>
  mutate(CollDate = as.Date(paste(month, '15', Year, sep='-'), format='%b-%d-%Y')) |>
  pivot_longer(cols=c('TN', `Inorganic N`, 'TP', 'Phosphate')) |>
  mutate(name = factor(name, levels=c('TN','TP','Inorganic N', 'Phosphate'))) |>
  mutate(WaterbodyName = factor(WaterbodyName, levels = c('Wind River Outlet', 'Muddy Creek', 'Fivemile Creek', 'Wind River Inlet')))

ggplot(trib_nuts) + 
  geom_point(aes(CollDate, value, color=WaterbodyName),size=2) +
  geom_line(aes(CollDate, value, color=WaterbodyName)) +
  scale_color_manual('',values =trib_colors) +
  facet_wrap(~name, scales = 'free') +
  labs(x='',y='Concentration '~(mg~L^-1)) +
  theme_minimal() +
  theme(legend.position = 'none')
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
  mutate(month=month(fakedate, label=TRUE)) |>
  group_by(year(fakedate), nutrient) |>
  mutate(sumpercRet=sum(perc_ret, na.rm=TRUE)) |>
  ungroup()


summary(lm(perc_ret~m, retention |> filter(nutrient=='TN')))
#slope=17.52
#p-value: 0.0007671
#Adjusted R-squared:  0.6267 

summary(lm(perc_ret~m, retention |> filter(nutrient=='TP')))
#slope=-8.924
#p-value: 0.4045
#Adjusted R-squared:  -0.00957 

summary(lm(perc_ret~m, retention |> filter(nutrient=='TP',
                                           perc_ret>-600)))
#slope=-17.110
#p-value: 0.01095
#Adjusted R-squared:0.1814 

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
  #geom_smooth(aes(m, perc_ret), method='reml', se=FALSE) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  facet_wrap(~nutrient, scales='free_y') +
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
  geom_point(aes(CollDate, concentration, fill=WaterbodyName), alpha=0.5, shape=21) +
  #geom_path(aes(CollDate, concentration, color=WaterbodyName), alpha=0.5) +
  scale_fill_viridis_d('', option='magma') +
  geom_point(aes(CollDate, meanConc),shape=22, size=2.5) +
  geom_line(aes(CollDate, meanConc, group=Year)) +
  facet_wrap(~nutrient, scales='free_y') +
  labs(x='',y='Nutrient concentration'~(mg~L^-1)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave('Figures/nutrientdynamics_boysen.png', height=4.5, width=6.5, units='in',dpi=1200)



# 9. Can we use discharge as a proxy for nutrient loading? ####

loadQ <- TribLoadFlux |>
  mutate(fakedate = paste0(Year, '-', month, '-01'),
         fakedate = as.Date(fakedate, format='%Y-%b-%d')) |>
  filter(!grepl('Outlet', WaterbodyName)) |>
  group_by(fakedate) |>
  summarise(discharge = sum(Discharge, na.rm=TRUE)) |>
  ungroup () |>
  left_join(TotLoad |>
              select(fakedate, nutrient, TotalLoad) |>
              filter(nutrient %in% c('TN','TP')) |>
              distinct()) |>
  # missing Sep 2023 discharge value from TribLoadFlux becuase there were no nutrient data that month. Add back in for use in PCA, NMDS later
  bind_rows(BoysenTribs |>
          filter(ShortName_Revised=='Discharge',
                 Year==2023,
                 !grepl('Outlet', WaterbodyName)) |>
          mutate(month=month(CollDate)) |>
          group_by(WaterbodyName, month, Year) |>
          summarise(Discharge=mean(ChemValue)) |>
          ungroup() |>
          filter(month==9) |>
          mutate(fakedate=as.Date('2023-09-01')) |>
          group_by(fakedate) |>
          summarise(discharge=sum(Discharge)) |>
          ungroup()) 



  

options(scipen = 999)
ggplot(loadQ, aes(discharge,TotalLoad, color =nutrient)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_minimal() +
  labs(x='Discharge'~(m^3~s^-1),
       y='Total nutrient load (kg)') +
  theme(legend.title=element_blank()) 

ggsave('Figures/dischargeNutrientLoad.png',width=6.5, height=4.5, units='in', dpi=1200)

n<-lm(TotalLoad~discharge, loadQ|>filter(nutrient=='TN'))
summary(n)

p<-lm(TotalLoad~discharge, loadQ|>filter(nutrient=='TP'))
summary(p)
# yes.








ggplot(TribLoadFlux |>
         mutate(month=factor(month, levels=c('Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')))) +
  geom_boxplot(aes(month, TP)) +

ggplot(TribLoadFlux |>
         mutate(month=factor(month, levels=c('Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')))) +
  geom_boxplot(aes(month, TN))










