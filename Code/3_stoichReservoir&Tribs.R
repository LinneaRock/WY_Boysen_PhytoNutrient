###########################
# Figure 3 nutrient ratios  
###########################

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')

# 2. get average N:P values ####
river_stoich <- TribLoadFlux |>
  mutate(TN.TP = (TN/TP)*2.11306) |> # makes molar from mg/L ratios
  select(WaterbodyName, Year, month, TN.TP) |>
  drop_na() |>
  mutate(type=ifelse(grepl('Outlet', WaterbodyName), 'Outlet', 'Tributary')) |>
  group_by(month, Year, type) |>
  mutate(meanNP=mean(TN.TP, na.rm=TRUE)) |>
  ungroup()


lake_stoich <- BoysenNutrient |>
  filter(ShortName_Revised=='TN.TP') |>
  pivot_wider(names_from = ShortName_Revised, values_from = ChemValue) |>
  select(WaterbodyName, Year, month, CollDate, TN.TP) |>
  drop_na() |>
  group_by(CollDate) |>
  mutate(meanNP=mean(TN.TP, na.rm=TRUE)) |>
  ungroup() |>
  mutate(type='Reservoir')

ave_stoich <- lake_stoich |>
  select(month, Year, meanNP) |>
  rename(Reservoir = meanNP) |>
  distinct() |>
  left_join(river_stoich |>
              select(month, Year, type, meanNP) |>
              distinct() |>
              pivot_wider(names_from = type, values_from = meanNP, id_cols = c(month, Year))) |>
  mutate(month = factor(month, levels=c('May','Jun','Jul','Aug','Sep','Oct')))


# 3. Plotting N:P comparisons  ####
ggplot(ave_stoich) +
  geom_point(aes(Tributary, Reservoir, color=month)) +
  geom_abline(slope=1,intercept=0) +
  scale_color_viridis_d('',option='magma') +
  theme_minimal() 


ggplot(ave_stoich) +
  geom_point(aes(Tributary, Outlet, color=month)) +
  geom_abline(slope=1,intercept=0) +
  scale_color_viridis_d('',option='magma') +
  theme_minimal()


ggplot(ave_stoich) +
  geom_point(aes(Reservoir, Outlet, color=month)) +
  geom_abline(slope=1,intercept=0) +
  scale_color_viridis_d('',option='magma') +
  theme_minimal()


boxplot_ave_df <- ave_stoich |>
  pivot_longer(cols=c(Reservoir, Outlet, Tributary), names_to = 'type', values_to = 'meanNP')  |>
  mutate(month = factor(month, levels=c('May','Jun','Jul','Aug','Sep','Oct'))) 


ggplot() +
  geom_boxplot(boxplot_ave_df, mapping=aes(type, meanNP,fill=month))  +
  scale_fill_viridis_d('',option='magma') +
  theme_minimal()


boxplot_df <- bind_rows(lake_stoich, river_stoich)  |>
  mutate(month = factor(month, levels=c('Jan', 'Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')),
         type = factor(type, levels=c('Tributary',"Reservoir",'Outlet'))) 

#~#
ggplot() +
  geom_boxplot(boxplot_df, 
               mapping=aes(type, TN.TP,fill=month))  +
  scale_fill_viridis_d('',option='magma') +
  theme_minimal() +
  labs(x='',y='TN:TP molar ratio')

ggsave('Figures/Fig3/fig3_NPboxplots.png',width=6.5,height=4.5,units='in',dpi=1200)
#~#  


ggplot() +
  geom_boxplot(boxplot_df |> 
                 filter(month %in% c('May','Jun','Jul','Aug','Sep','Oct')), 
               mapping=aes(type, TN.TP,fill=month))  +
  scale_fill_viridis_d('',option='magma') +
  theme_minimal()

ggplot() +
  geom_boxplot(boxplot_df |> 
                 filter(!month %in% c('May','Jun','Jul','Aug','Sep','Oct')), 
               mapping=aes(type, TN.TP,fill=month))  +
  scale_fill_viridis_d('',option='magma') +
  theme_minimal()

