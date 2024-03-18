#-----------------------------------------------#
# copying Kelsey's code
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. Prep data for plotting ####
BoysenChemPhys <- ChemPhys |>
  dplyr::select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  mutate(ChemUnits = sub('/l', '/L', ChemUnits)) |>
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  dplyr::select(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, BelowDet, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct() |>
  mutate(year = year(CollDate)) |>
  filter(year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) 




BoysenPhyto <- BoysenChemPhys |>
  dplyr::select(StationID, WaterbodyName)  |>
  left_join(Phyto|>
              rename(chemName = WaterbodyName)) |>
  drop_na(CollDate) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  unique() |>
  group_by(StationID, WaterbodyName, CollDate, month, Year, `Genus/Species/Variety`) |>
  summarise(indsum = sum(`Individuals (Raw Cnt)`)) |>
  ungroup() |>
  group_by(StationID, WaterbodyName, month, Year) |>
  mutate(totaln=sum(indsum)) |>
  ungroup() |>
  mutate(perc =(indsum/totaln)*100) |>
  group_by(StationID, month, Year) |>
  mutate(check = sum(perc)) |>
  ungroup()

ggplot(BoysenPhyto) +
  geom_bar(aes(month, perc, fill=`Genus/Species/Variety`), stat='identity') +
  facet_wrap(WaterbodyName~Year)

n_phyto <- BoysenPhyto |>
  select(`Genus/Species/Variety`) |>
  distinct()


# 2. BETA DIVERSITY ANALYSIS ####
#PERMANOVA with adonis or adonis2
#can help test relationships between water quality values (for example) and community composition 
library(vegan)


# I need all my other variables with a Grouping variable column to match with the dist data
sd_data <- BoysenChemPhys |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN',
                                       ShortName_Revised=='Total Ammonia as N'~'NH4',
                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                       ShortName_Revised=='Orthophosphate as P (total)'~'PO4',
                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                       ShortName_Revised=='Chlorophyll a (phytoplankton)'~'CHLA')) |>
  select(Group, WaterbodyName, CollDate, Latitude, Longitude, ShortName_Revised, ChemValue) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue)

phyto_data <- BoysenPhyto |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, `Genus/Species/Variety`, indsum, WaterbodyName, CollDate, Year, month) |>
  left_join(sd_data) |>
  pivot_wider(names_from = `Genus/Species/Variety`, values_from = indsum) |>
  as.data.frame() 

# create metadata 
sd_phyto <- phyto_data |>
  select(Group, WaterbodyName, CollDate, Year, month, Latitude, Longitude, PO4, NH4,TP, TN, NO3, CHLA)

# create distance matrix
dist_phyto <- phyto_data |>
  select(-WaterbodyName, -CollDate, -Latitude, -Longitude, -PO4, -NH4,-TP, -TN, -NO3, -CHLA, -Year, -month)
rownames(dist_phyto) <- dist_phyto$Group
dist_phyto <- dist_phyto[,-1]
dist_phyto <- as.matrix(dist_phyto)
dist_phyto <- replace(dist_phyto, is.na(dist_phyto), 0)

dist <- vegdist(dist_phyto, method = 'bray')

dist

adonis2(dist~Latitude, sd_phyto)
adonis2(dist~WaterbodyName, sd_phyto, strata=sd_phyto$CollDate) #??
adonis2(dist~CHLA, sd_phyto) # chla not collected 2021-05-18
adonis2(dist~TN, sd_phyto)
adonis2(dist~TP, sd_phyto)
adonis2(dist~TN + TP, sd_phyto) # ??
adonis2(dist~NO3, sd_phyto)
adonis2(dist~NH4, sd_phyto)
adonis2(dist~PO4, sd_phyto)



# 3. NMDS ####
nmds <- metaMDS(dist)

nmds # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good

scores <- scores(nmds) |>
  as_tibble(rownames='Group') |>
  left_join(sd_phyto)


ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=WaterbodyName))

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=as.character(Year)))

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=CHLA)) +
  scale_color_viridis_c()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=TN)) +
  scale_color_viridis_c()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=TP)) +
  scale_color_viridis_c()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=PO4)) +
  scale_color_viridis_c()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=NO3)) +
  scale_color_viridis_c()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=NH4)) +
  scale_color_viridis_c()



ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=month, shape=as.character(Year))) +
  geom_line(aes(group=WaterbodyName), alpha=0.25)

#You can also investigate the species which may be driving the site distribution pattern, referred to as intrinsic variables.
spp.fit <- envfit(nmds, dist_phyto, permutations=999)
head(spp.fit)

# from Jordy 
# check out ordiplot()
# i can use lines to connect my points in nmds by site and time


ordiplot(nmds, type='n', main='intrinsic species')
plot(spp.fit, p.max = 0.01, col = "black", cex = 0.7) # change the significance level of species shown with p.max

# Environmental variables can also be used with envfit which are referred to as extrinsic variables. This works best with continuous variables of which there is only one (A1) in this dataset.If you only want to fit vector variables (continuous variables) use vectorfit and if you only want to fit factor variables (categorical variables) use factorfit but envfit can do this automatically.

env.fit <- envfit(nmds, sd_phyto, permutations=999, na.rm=TRUE)
head(env.fit)

ordiplot(nmds, type='n', main='extrinsic factors')
plot(env.fit, p.max = 0.01, col = "black", cex = 0.7) 

# 4. BETA DIVERSITY DISPERSION ####
#beta Dispersion plot

#transformation.. a critical 1st step before looking at beta diversity.. Need to account for differences in the number of counts in each sample before doing any distance matrices. 











































