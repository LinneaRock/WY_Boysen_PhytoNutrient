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
         Spp = `Genus/Species/Variety`) |>
  filter(CommonName == 'Boysen Reservoir')
#table(combine$Genus.Species.Variety)

# count model 
m1 <- glm(data = combine, Count ~ TN_mgL + Elevation, family = 'poisson')
summary(m1)

library(ggeffects)
newd <- data.frame(ggpredict(m1, terms = c('TN_mgL[all]', 'Elevation[all]'), interval = 'confidence', ci.lvl = 0.9))

ggplot(newd, aes(x, predicted)) +
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = facet), alpha = 0.6,
  #             linetype = 'dashed', color = 'black', linewidth = 1.2) +
  geom_boxplot(aes(color = group))


