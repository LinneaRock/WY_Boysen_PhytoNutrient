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
  mutate(Year = year(CollDate)) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) 

phyto_class <- read.csv('Data/phyto_class.csv') |>
  rename(cat=Class) 

phyto_class$cat <- sub('^$', 'unknown', phyto_class$cat)


BoysenPhytoA <- BoysenChemPhys |>
  dplyr::select(StationID, WaterbodyName)  |>
  left_join(Phyto|>
              rename(chemName = WaterbodyName,
            Genus.Species.Variety = `Genus/Species/Variety`)) |>
  drop_na(CollDate) |>
  filter(Year >= 2020) |>
  mutate(month = month(CollDate, label=TRUE, abbr=TRUE)) |>
  unique() |>
  left_join(phyto_class) 
  



BoysenPhyto <- BoysenPhytoA |>
  filter(RepNum == 0) |> #only keep first counts  to avoid the insane confusion I had when I ignored this column
  distinct()|>
  group_by(WaterbodyName, CollDate, month, Year, Genus.Species.Variety) |>
  mutate(indsum = sum(`Individuals (Raw Cnt)`)) |>
  ungroup() |>
  distinct() |>
  group_by(WaterbodyName, month, Year) |>
  mutate(totaln=sum(indsum)) |>
  ungroup() |>
  distinct() |>
  mutate(perc =(indsum/totaln)*100) |>
  group_by(WaterbodyName, month, Year) |>
  mutate(check = sum(perc)) |>
  ungroup() |>
  select(WaterbodyName, CollDate, month, Year, Genus.Species.Variety, indsum, totaln, perc, check)



BoysenPhyto_cat <- BoysenPhytoA |>
  filter(RepNum == 0) |>
  left_join(BoysenPhyto |> select(WaterbodyName, Year, month, totaln)) |>
  distinct() |>
  group_by(WaterbodyName, CollDate, month, Year, cat, totaln) |>
  summarise(catsum = sum(`Individuals (Raw Cnt)`)) |>
  ungroup() |>
  mutate(Catperc =(catsum/totaln)*100) |>
  group_by(WaterbodyName, month, Year) |>
  mutate(checkCat = sum(Catperc)) |>
  ungroup() |>
  select(WaterbodyName, month, Year, cat, Catperc) |>
  distinct() |>
  pivot_wider(names_from = cat, values_from = Catperc) 
BoysenPhyto_cat <- replace(BoysenPhyto_cat, is.na(BoysenPhyto_cat), 0)





ggplot(BoysenPhyto) +
  geom_bar(aes(month, perc, fill=Genus.Species.Variety), stat='identity') +
  facet_wrap(WaterbodyName~Year)

n_phyto <- BoysenPhyto |>
  select(Genus.Species.Variety) |>
  distinct()


# 2. BETA DIVERSITY ANALYSIS ####
#PERMANOVA with adonis or adonis2
#can help test relationships between water quality values (for example) and community composition 
library(vegan)


# I need all my other variables with a Grouping variable column to match with the dist data
# create metadata 
sd_data <- BoysenChemPhys |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN',
                                       ShortName_Revised=='Total Ammonia as N'~'NH4',
                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                       ShortName_Revised=='Orthophosphate as P (total)'~'PO4',
                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                       ShortName_Revised=='Chlorophyll a (phytoplankton)'~'CHLA')) |>
  select(Group, WaterbodyName, CollDate, Year, month, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`, unknown, Crustacean, Cyanobacteria, Dinoflagellate) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue)


# match up data for later
phyto_data <- BoysenPhyto |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, Genus.Species.Variety, indsum, WaterbodyName, CollDate, Year, month) |>
  left_join(sd_data) |>
  pivot_wider(names_from = Genus.Species.Variety, values_from = indsum) |>
  as.data.frame() 


# create distance matrix
dist_phyto <- phyto_data |>
  select(-WaterbodyName, -CollDate, -Latitude, -Longitude, -PO4, -NH4,-TP, -TN, -NO3, -CHLA, -Year, -month, -Diatom, -`Green algae`, -unknown, -Crustacean, -Cyanobacteria, -Dinoflagellate)
rownames(dist_phyto) <- dist_phyto$Group
dist_phyto <- dist_phyto[,-1]
dist_phyto <- as.matrix(dist_phyto)
dist_phyto <- replace(dist_phyto, is.na(dist_phyto), 0)

dist <- vegdist(dist_phyto, method = 'bray')

dist

adonis2(dist~Latitude, sd_data)
adonis2(dist~WaterbodyName, sd_data, strata=sd_data$CollDate) #??
adonis2(dist~WaterbodyName, sd_data) # high p-value == sites are the same in terms of their beta diversity (i.e., comparing samples to each other and answers question 'how different')? 
adonis2(dist~CHLA, sd_data) # chla not collected 2021-05-18
adonis2(dist~TN, sd_data)
adonis2(dist~TP, sd_data)
adonis2(dist~TN * TP, sd_data) 
adonis2(dist~NO3, sd_data)
adonis2(dist~NH4, sd_data) # produces significant result
adonis2(dist~PO4, sd_data) # produces significant result
adonis2(dist~Cyanobacteria, sd_data) # produces significant result



# 3. NMDS ####
nmds <- metaMDS(dist)

nmds # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good

scores <- scores(nmds) |>
  as_tibble(rownames='Group') |>
  left_join(sd_data)


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
  geom_point(aes(color=Cyanobacteria)) +
  scale_color_viridis_c()

# from Jordy 
# check out ordiplot()
# i can use lines to connect my points in nmds by site and time

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=month, shape=as.character(Year))) +
  geom_line(aes(group=WaterbodyName), alpha=0.25)

#You can also investigate the species which may be driving the site distribution pattern, referred to as intrinsic variables.
spp.fit <- envfit(nmds, dist_phyto, permutations=999)
head(spp.fit)


# ordiplot(nmds, type='n', main='intrinsic species')
# plot(spp.fit, p.max = 0.01, col = "black", cex = 0.7) # change the significance level of species shown with p.max


spp.fit_df <- as.data.frame(scores(spp.fit, display='vectors'))  #extracts relevant scores from sppifit
spp.fit_df <- cbind(spp.fit_df, spp.variables = rownames(spp.fit_df)) #and then gives them their names

spp.fit_df <- cbind(spp.fit_df, pval = spp.fit$vectors$pvals) # add pvalues to dataframe
sig.spp.fit <- subset(spp.fit_df, pval<=0.01) #subset data to show variables significant at 0.05

sig.spp.fit

#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!

ggplot() +
  geom_point(scores, mapping=aes(x=NMDS1, y=NMDS2, color=Cyanobacteria)) +
  theme_minimal() +
  scale_color_viridis_c() +
  geom_segment(sig.spp.fit, mapping=aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(sig.spp.fit, mapping=aes(x=NMDS1, y=NMDS2, label = spp.variables), cex = 3, direction = "both", segment.size = 0.25) #add labels, use ggrepel::geom_text_repel so that labels do not overlap






# Environmental variables can also be used with envfit which are referred to as extrinsic variables. This works best with continuous variables of which there is only one (A1) in this dataset.If you only want to fit vector variables (continuous variables) use vectorfit and if you only want to fit factor variables (categorical variables) use factorfit but envfit can do this automatically.

env.fit <- envfit(nmds, sd_data, permutations=999, na.rm=TRUE)
head(env.fit)
# ordiplot(nmds, type='n', main='extrinsic factors')
# plot(env.fit, p.max = 0.01, col = "black", cex = 0.7) 



env.fit_df <- as.data.frame(scores(env.fit, display='vectors'))  #extracts relevant scores from envifit
env.fit_df <- cbind(env.fit_df, env.variables = rownames(env.fit_df)) #and then gives them their names

env.fit_df <- cbind(env.fit_df, pval = env.fit$vectors$pvals) # add pvalues to dataframe
sig.env.fit <- subset(env.fit_df, pval<=0.05) #subset data to show variables significant at 0.05

sig.env.fit

#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!

ggplot() +
  geom_point(scores, mapping=aes(x=NMDS1, y=NMDS2, color=Cyanobacteria)) +
  theme_minimal() +
  scale_color_viridis_c() +
  geom_segment(sig.env.fit, mapping=aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(sig.env.fit, mapping=aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 3, direction = "both", segment.size = 0.25) #add labels, use ggrepel::geom_text_repel so that labels do not overlap



# 4. BETA DIVERSITY DISPERSION ####
#beta Dispersion plot --- I think these are the NMDS plots??????



# 5. #ALPHA DIVERSITY ####
#Taxonomy Line Plot
#rarefy and normalizing before calculating any alpha diversity metrics- since my samples have a variety of number of sequences or organisms


# 7. %cyano plot from Smith 1983 ####

# multiply N:P mass by 2.11306 to get molar
NP_cyano <- sd_data |>
  mutate(TN.TP = (TN/TP)*2.11306,
         NO3.PO4 = (NO3/PO4)*2.11306,
         NH4.PO4 = (NH4/PO4)*2.11306,
         IN.PO4 = ((NH4+NO3)/PO4)*2.11306) |>
  pivot_longer(cols=c(TN.TP, NO3.PO4, NH4.PO4, IN.PO4), names_to = 'NP_type')

ggplot(NP_cyano) +
  geom_point(aes(value, Cyanobacteria, color=NP_type)) +
  labs(x='N:P molar ratio', y='% Cyanobacteria abundance') +
  geom_vline(xintercept = 29) +
  theme_minimal()
  

# 8. Tributary Budget Plots ####
bud <- BoysenTribs |>
  select(-X) |>
  mutate(month=month(CollDate, label=TRUE, abbr=TRUE),
         Year=year(CollDate)) |>
  filter(Year>=2020) |>
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN',
                                       ShortName_Revised=='Total Ammonia as N'~'NH4',
                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                       ShortName_Revised=='Dissolved Orthophosphate'~'PO4',
                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                       ShortName_Revised=='Chlorophyll a (phytoplankton)'~'CHLA',
                                       ShortName_Revised==''))
  pivot_wider(idcols=c(WaterbodyName, CollDate),names_from = ShortName_Revised, values_from = ChemValue)
  







































