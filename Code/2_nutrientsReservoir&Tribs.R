##########################################################
# Figure 2 nutrient loading and reservoir concentrations 
##########################################################

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')

# dissolved inorganic nitrate
TribLoadFlux |> 
  mutate(IN_load_kg = NH4_load_kg + NO3_load_kg) |>
  mutate(TotalTrib_IN_kg = TotalTrib_NH4_kg + TotalTrib_NO3_kg) |>
  mutate(IN_load_kg=ifelse(WaterbodyName=='Wind River Outlet', -1*IN_load_kg, NA)) |>
  mutate(fakedate = as.Date(paste0(Year,'-',month,'-01'), format='%Y-%b-%d')) |>
  select(fakedate, TotalTrib_IN_kg, IN_load_kg) |>
  distinct() |>
  mutate(dummy = ifelse(is.na(IN_load_kg), 'a','b')) |>
  arrange(dummy,fakedate) |>
  filter(!if_all(c(TotalTrib_IN_kg, IN_load_kg), is.na)) |>
  ggplot() +
  geom_point(aes(fakedate, TotalTrib_IN_kg),size=3) +
  geom_path(aes(fakedate, TotalTrib_IN_kg)) +
  geom_point(aes(fakedate, IN_load_kg),size=3) +
  geom_path(aes(fakedate, IN_load_kg)) +
  geom_abline(slope=0, intercept=0) +
  theme_minimal() +
  annotate("rect", xmin = as.Date('2020-01-01'), xmax = as.Date('2020-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
  annotate("rect", xmin = as.Date('2020-10-15'), xmax = as.Date('2021-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
  annotate("rect", xmin = as.Date('2021-10-15'), xmax = as.Date('2022-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
  annotate("rect", xmin = as.Date('2022-10-15'), xmax = as.Date('2023-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
  labs(x='',y='Nitrate + ammonium load (kg)')
ggsave('Figures/ASLO24/totalIN_loading.png',height = 4.5, width = 6.5, units='in', dpi=1200)
