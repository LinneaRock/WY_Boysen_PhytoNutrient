#---------------------------#
# Boysen Profile Comparisons
#---------------------------#

source('Data/CALL_DATA_LIB.R')
library(vegan)


#1. create distance matrix ####
dist_temp<- BoysenProfile  |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, depth_m, temp_C) |>
  pivot_wider(id_cols = Group, names_from = depth_m, values_from = temp_C)
  
  
  
  
  
  select(-WaterbodyName, -CollDate, -julianday, -Latitude, -Longitude, -PO4, -NH4,-TP, -TN, -NO3, -CHLA, -TN.TP,  -IN.PO4, -Year, -month, -Diatom, -`Green algae`, -Flagellate, -`Golden algae`, -Cyanobacteria, -Dinoflagellate, -DO, -Secchi, -pH, -SpC,-H)
rownames(dist_phyto) <- dist_phyto$Group
dist_phyto <- dist_phyto[,-1]
dist_phyto <- as.matrix(dist_phyto)
dist_phyto <- replace(dist_phyto, is.na(dist_phyto), 0)


dist <- vegdist(dist_phyto, method = 'bray')

dist
