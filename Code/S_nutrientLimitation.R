#-------------------------------------------------------------------#
# Assessing nutrient limitation via some of my favorite est. methods 
#-------------------------------------------------------------------#

source('Data/CALL_DATA_LIB.R')

lim_dat <- BoysenNutrient |>
  select(-ChemUnits, -julianday, -SampleDepth) |>
  pivot_wider(names_from = ShortName_Revised, values_from = ChemValue) |>
  mutate(DIN.TP = ((NO3+NH4)/TP) * 2.11306)  |> # multiply by constant to make molar ratio
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) 


# TN:TP limitation ####
ggplot() +
  geom_jitter(lim_dat, mapping = aes(log10(TP), log10(TN.TP), 
                                    color = WaterbodyName), size=2) +
  geom_abline(slope = 0, intercept = log10(41), linetype = "dashed") + # bergstrom P limitation line
  geom_abline(slope = 0, intercept = log10(19), linetype = "dashed") +  # bergstrom N limitation line
  geom_abline(slope = 0, intercept = log(16, base = 10), color = "red4") + # redfield
  geom_vline(xintercept = log(30, base = 10)) + # dodds mccauley N limitation P > 30
  geom_abline(slope = 0, intercept = log(32, base = 10)) + # dodds mccauley N limitation TN:TP < 14
  geom_abline(slope = 0, intercept = log(38, base = 10), color = '#336a98') + # Sakamoto, 1966; Smith 1982; Rhee 1980, Forsberg 1980  
  geom_abline(slope = 0, intercept = log(22, base = 10), color = '#336a98') + # Sakamoto, 1966; Smith 1982; Rhee 1980, Forsberg 1980
  geom_abline(slope = 0, intercept = log(53, base = 10), color = "#ffc857") + # Ptacnik, 2010
  theme_minimal() +
  labs(y = "Log TN:TP", x = "Log TP") +
  annotate('text', label = 'Redfield 16:1 line', x = -2.5, y = 1.15, hjust = 0, size = 2.5, color = "red4") +
  annotate('text', label = 'Predicted N limitation\n below dashed line \n (Bergström, 2010)', 
           x = 0.55, y = 1, hjust = 0, size = 2.5) +
  annotate('text', label = 'Predicted P limitation\n above dashed line \n (Bergström, 2010)', 
           x = 0.5, y = 1.85, hjust = 0, size = 2.5) +
  annotate('text', label = 'Predicted N limitation below blue line \n (Forsberg, 1980; Rhee, 1980; \n Sakamoto, 1966; Smith, 1982)', x = -0.9, y = 1.1, hjust = 0, size = 2.5, color = '#336a98') +
  annotate('text', label = 'Predicted P limitation above blue line \n (Forsberg, 1980; Rhee, 1980; \n Sakamoto, 1966; Smith, 1982)', x = -1.8, y = 1.85, hjust = 0, size = 2.5, color = '#336a98') + 
  annotate('text', label = 'Predicted N limitation \n left of black vertical line, \n below black horizontal line \n (Downing & McCauley, 1992)', 
           x = -0.5, y = 0.75, hjust = 0, size = 2.5) +
  annotate('text', label = 'Predicted P limitation \n (Ptacnik et al., 2010)',
           x = 0.55, y = 2.05, hjust = 0, size = 2.5, color = "#ffc857") +
  scale_color_viridis_d('',option='magma') +
  theme(legend.position = 'none')

ggsave('Figures/nutrient_limitation.png', height=5,width=7,units='in',dpi=1200)










# TN:TP limitation by month ####
ggplot() +
  geom_jitter(lim_dat, mapping = aes(TP, TN.TP, 
                                     color = month), size=2) +
  geom_abline(slope = 0, intercept = 41, linetype = "dashed") + # bergstrom P limitation line
  geom_abline(slope = 0, intercept = 19, linetype = "dashed") +  # bergstrom N limitation line
  geom_abline(slope = 0, intercept = 16, color = "red4") + # redfield
 # geom_vline(xintercept = 30) + # dodds mccauley N limitation P > 30
  geom_abline(slope = 0, intercept = 32) + # dodds mccauley N limitation TN:TP < 14
  geom_abline(slope = 0, intercept = 38, color = '#336a98') + # Sakamoto, 1966; Smith 1982; Rhee 1980, Forsberg 1980  
  geom_abline(slope = 0, intercept = 22, color = '#336a98') + # Sakamoto, 1966; Smith 1982; Rhee 1980, Forsberg 1980
  geom_abline(slope = 0, intercept = 53, color = "#ffc857") + # Ptacnik, 2010
  theme_minimal()# +
  # labs(y = "Log TN:TP", x = "Log TP") +
  # annotate('text', label = 'Redfield 16:1 line', x = -2.5, y = 1.15, hjust = 0, size = 2.5, color = "red4") +
  # annotate('text', label = 'Predicted N limitation\n below dashed line \n (Bergström, 2010)', 
  #          x = 0.55, y = 1, hjust = 0, size = 2.5) +
  # annotate('text', label = 'Predicted P limitation\n above dashed line \n (Bergström, 2010)', 
  #          x = 0.5, y = 1.85, hjust = 0, size = 2.5) +
  # annotate('text', label = 'Predicted N limitation below blue line \n (Forsberg, 1980; Rhee, 1980; \n Sakamoto, 1966; Smith, 1982)', x = -0.9, y = 1.1, hjust = 0, size = 2.5, color = '#336a98') +
  # annotate('text', label = 'Predicted P limitation above blue line \n (Forsberg, 1980; Rhee, 1980; \n Sakamoto, 1966; Smith, 1982)', x = -1.8, y = 1.85, hjust = 0, size = 2.5, color = '#336a98') + 
  # annotate('text', label = 'Predicted N limitation \n left of black vertical line, \n below black horizontal line \n (Dodds & McCauley, 1992)', 
  #          x = -0.5, y = 0.75, hjust = 0, size = 2.5) +
  # annotate('text', label = 'Predicted P limitation \n (Ptacnik et al., 2010)',
  #          x = 0.55, y = 2.05, hjust = 0, size = 2.5, color = "#ffc857") +
  # scale_color_viridis_d('',option='turbo') 








# DIN:TP limitation ####
# maybe not viable metric because of DIN detection limits
ggplot() +
  geom_jitter(lim_dat, mapping = aes(log10(TP), log10(DIN.TP), 
                                     color = WaterbodyName), size=2) +
  geom_abline(slope = 0, intercept = log10(3.4), linetype = "dashed") + # bergstrom P limitation line
  geom_abline(slope = 0, intercept = log10(1.5), linetype = "dashed") +  # bergstrom N limitation line
  theme_minimal() +
  labs(y = "Log DIN:TP", x = "Log TP") +
  annotate('text', label = 'Predicted N limitation below dashed line \n (Bergström, 2010)', 
           x = -2.5, y = 0.17, hjust = 0, size = 4) +
  annotate('text', label = 'Predicted P limitation above dashed line \n (Bergström, 2010)', 
           x = -2.5, y = 0.55, hjust = 0, size = 4) +
  scale_color_viridis_d('',option='turbo')

