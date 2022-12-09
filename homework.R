# looking at phyto counts for DataViz homework Dec. 2022

library(tidyverse)

phyto.count <- read_csv('Data/RawData_WYDEQ/Phytoplankton_2013-2021.csv', skip = 5)
chemphys <- read_csv("Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv")

phyto.count2 <- phyto.count |>
  select('StationID', 'CollDate', 'Genus/Species/Variety', 'Individuals (Raw Cnt)') |>
  distinct() 

chemphys2 <- chemphys |>
  filter(SampleDepth == 0.5,
         ShortName_Revised %in% c('Phosphorus as P (total)', 'Total Nitrogen (unfiltered)')) |>
  pivot_wider(names_from = c('ShortName_Revised', 'ChemUnits'), values_from = 'ChemValue') |>
  distinct() 
  

combine <- inner_join(phyto.count2, chemphys2,  by = c("StationID", "CollDate")) |>
  rename(TN_mgL = `Total Nitrogen (unfiltered)_mg/l`,
         TP_mgL = `Phosphorus as P (total)_mg/l`,
         Count = `Individuals (Raw Cnt)`,
         Spp = `Genus/Species/Variety`) #|>
  #filter(CommonName == 'Boysen Reservoir')
#table(combine$Genus.Species.Variety)

# count model 
m1 <- glm(data = combine, Count ~ TN_mgL + Spp + Elevation, family = 'poisson')
summary(m1)

library(ggeffects)
# how do counts of a single zoop species vary with increasing TN concentration and how does this relationship vary across elevation?
newd <- data.frame(ggpredict(m1, terms = c('TN_mgL[all]', 'Spp[Dictyosphaerium]', 'Elevation[all]'), interval = 'confidence', ci.lvl = 0.9)) |>
  mutate(facet = as.factor(facet))

ggplot(newd, aes(x, predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = facet), alpha = 0.6,
                linetype = 'dashed', color = 'black', linewidth = 1.2) +
  geom_line(aes(color = facet, linetype = facet), linewidth = 1.2) +
  scale_fill_viridis_d(option = 'viridis', begin = 0.2, end = 0.5, direction = 1, # These two scale_fill/colors with all the same specs provides one nice legend!!
                       name = '',
                       breaks = c('4726.48','9050'),
                       labels = c('Plains Lake','Alpine Lake')) +
  scale_color_viridis_d(option = 'viridis', begin = 0.2, end = 0.5, direction = 1,
                        name = '',
                        labels = c('Plains Lake', 'Alpine Lake'),
                        breaks = c('4726.48','9050')) +
  labs(x = 'Total Nitrogen'~mg~L^-1, y = "Dictyosphaerium count") +
  guides(linetype = 'none') +
  theme_classic(base_size = 12) +
  theme(text = element_text(family = 'serif'),
        axis.text = element_text(color = 'black'))

ggsave('zoopcounthomework.png', height = 4.5, width = 6.5, units = 'in', dpi = 1200)


# alternative count model 
m1 <- glm(data = combine, Count ~ TN_mgL + Spp + Elevation, family = 'poisson')
summary(m1)

library(ggeffects)
# how do counts of zoop species vary with increasing TN concentration and how does this relationship vary across elevation?
newd <- data.frame(ggpredict(m1, terms = c('TN_mgL[all]', 'Spp[Dictyosphaerium,Aphanizomenon,Hippodonta sp.,Oocystis]', 'Elevation[all]'), interval = 'confidence', ci.lvl = 0.9)) |>
  mutate(facet = as.factor(facet))

ggplot(newd, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = facet), alpha = 0.6,
              linetype = 'dashed', color = 'black', linewidth = 1.2) +
  geom_line(aes(color = facet, linetype = facet), linewidth = 1.2)  +
  facet_wrap(.~group, scales = 'free_y') +
  scale_fill_viridis_d(option = 'viridis', begin = 0.2, end = 0.5, direction = 1, # These two scale_fill/colors with all the same specs provides one nice legend!!
                       name = '',
                       breaks = c('4726.48','9050'),
                       labels = c('Plains Lake','Alpine Lake')) +
  scale_color_viridis_d(option = 'viridis', begin = 0.2, end = 0.5, direction = 1,
                        name = '',
                        labels = c('Plains Lake', 'Alpine Lake'),
                        breaks = c('4726.48','9050')) +
  labs(x = 'Total Nitrogen'~mg~L^-1, y = "Dictyosphaerium count") +
  guides(linetype = 'none') +
  theme_classic(base_size = 12) +
  theme(text = element_text(family = 'serif'),
        axis.text = element_text(color = 'black'))

ggsave('zoopcounthomework_option2.png', height = 4.5, width = 6.5, units = 'in', dpi = 1200)
