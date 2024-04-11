#-----------------------------------------------#
# copying Kelsey's code and Jordy's code and doing other things
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. Prep data for plotting ####
head(BoysenNutrient)
head(BoysenPhyto) 
head(BoysenTribs)

phyto_class <- read.csv('Data/phyto_class.csv') |>
  rename(cat=Class) 

phyto_class$cat <- sub('^$', 'unknown', phyto_class$cat)



BoysenPhyto_A <- BoysenPhyto |>
  # get rid of the zoops
  filter(Division != 'Copepoda',
         !Genus.Species.Variety %in% c('Daphnia','Cladocera', '	
Copepoda: nauplius')) |>
  left_join(phyto_class)  |>
  filter(RepNum == 0) |> #only keep first counts  to avoid the insane confusion I had when I ignored this column
  distinct()|>
  group_by(WaterbodyName, CollDate, month, Year, Genus.Species.Variety) |>
  mutate(indsum = sum(`Individuals (Raw Cnt)`),
         biovolumeSum_cellsL=sum(`Density (cells/L)`)) |>
  ungroup() |>
  distinct() |>
  group_by(WaterbodyName, month, Year) |>
  mutate(totaln=sum(indsum),
         totalbiovolume_cellsL=sum(biovolumeSum_cellsL)) |>
  ungroup() |>
  distinct() |>
  mutate(percCount =(indsum/totaln)*100,
         percBiovolume=(biovolumeSum_cellsL/totalbiovolume_cellsL)*100) |>
  group_by(WaterbodyName, month, Year) |>
  mutate(checkn = sum(percCount),
         checkbv=sum(percBiovolume)) |>
  ungroup() |>
  select(WaterbodyName, CollDate, month, Year, Genus.Species.Variety, indsum, totaln, percCount, checkn,biovolumeSum_cellsL, totalbiovolume_cellsL,percBiovolume,checkbv)
# percent biovolume = percent count -- duh but needed to verify



BoysenPhyto_cat <- BoysenPhyto |>
  # get rid of the zoops
  filter(Division != 'Copepoda',
         !Genus.Species.Variety %in% c('Daphnia','Cladocera', '	
Copepoda: nauplius')) |>
  filter(RepNum == 0) |>
  left_join(phyto_class)  |>
  left_join(BoysenPhyto_A |> select(WaterbodyName, Year, month, totaln, totalbiovolume_cellsL)) |>
  distinct() |>
  group_by(WaterbodyName, CollDate, month, Year, cat, totaln, totalbiovolume_cellsL) |>
  summarise(catsum = sum(`Individuals (Raw Cnt)`),
            CATbiovolumeSum_cellsL=sum(`Density (cells/L)`)) |>
  ungroup() |>
  mutate(Catperc =(catsum/totaln)*100,
         CATpercBiovolume=(CATbiovolumeSum_cellsL/totalbiovolume_cellsL)*100) |>
  group_by(WaterbodyName, month, Year) |>
  mutate(checkCat = sum(Catperc),
         checkbv=sum(CATpercBiovolume)) |>
  # percent biovolume = percent count -- duh but needed to verify
  ungroup() |>
  select(WaterbodyName, month, Year, cat, Catperc) |>
  distinct() |>
  pivot_wider(names_from = cat, values_from = Catperc) 
BoysenPhyto_cat <- replace(BoysenPhyto_cat, is.na(BoysenPhyto_cat), 0)



ggplot(BoysenPhyto_A) +
  geom_bar(aes(month, percCount, fill=Genus.Species.Variety), stat='identity') +
  facet_wrap(WaterbodyName~Year)

n_phyto <- BoysenPhyto_A |>
  select(Genus.Species.Variety) |>
  distinct() # 71 species found


# 2. BETA DIVERSITY ANALYSIS ####
#PERMANOVA with adonis or adonis2
#can help test relationships between water quality values (for example) and community composition 
library(vegan)


# I need all my other variables with a Grouping variable column to match with the dist data
# create metadata 
sd_data <- BoysenNutrient |>
  rbind(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN',
                                       ShortName_Revised=='Total Ammonia as N'~'NH4',
                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                       ShortName_Revised=='Orthophosphate as P (total)'~'PO4',
                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                       ShortName_Revised=='Chlorophyll a (phytoplankton)'~'CHLA',
                                       ShortName_Revised == 'Secchi Depth' ~'Secchi',
                                       ShortName_Revised=='pH'~'pH',
                                       ShortName_Revised=='DO, mg/L'~'DO')) |>
  filter(!is.na(ShortName_Revised)) |>
  select(Group, WaterbodyName, CollDate, Year, month, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) |>
  mutate(TN.TP = (TN/TP)*2.11306,
         NO3.PO4 = (NO3/PO4)*2.11306,
         NH4.PO4 = (NH4/PO4)*2.11306,
         IN.PO4 = ((NH4+NO3)/PO4)*2.11306)


# match up data for later
phyto_data <- BoysenPhyto_A |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, Genus.Species.Variety, indsum, WaterbodyName, CollDate, Year, month) |>
  left_join(sd_data) |>
  pivot_wider(names_from = Genus.Species.Variety, values_from = indsum) |>
  as.data.frame() 


# create distance matrix
dist_phyto <- phyto_data |>
  select(-WaterbodyName, -CollDate, -Latitude, -Longitude, -PO4, -NH4,-TP, -TN, -NO3, -CHLA, -TN.TP, -NO3.PO4, -NH4.PO4, -IN.PO4, -Year, -month, -Diatom, -`Green algae`, -Flagellate, -`Golden algae`, -Cyanobacteria, -Dinoflagellate, -DO, -Secchi, -pH)
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
adonis2(dist~TN*TP, sd_data)
adonis2(dist~TN.TP, sd_data)
adonis2(dist~NO3, sd_data)
adonis2(dist~NH4, sd_data) # produces significant result
adonis2(dist~PO4, sd_data) # produces significant result
adonis2(dist~NO3.PO4, sd_data) # sig
adonis2(dist~NH4.PO4, sd_data) # sig
adonis2(dist~IN.PO4, sd_data) # sig
adonis2(dist~Cyanobacteria, sd_data) # produces significant result
#adonis2(dist~TN.TP+NO3.PO4+NH4.PO4+IN.PO4+TN*TP*NO3*NH4*PO4, sd_data) # interesting result potentially here with interactions of N and P



# 3. NMDS and 4. BETA DIVERSITY DISPERSION ####
#beta Dispersion plot --- I think these are the NMDS plots?????? 
nmds <- metaMDS(dist)

nmds # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good

scores <- scores(nmds) |>
  as_tibble(rownames='Group') |>
  left_join(sd_data)


ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=WaterbodyName)) +
  theme_minimal()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=month)) +
  theme_minimal()

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
  scale_color_viridis_c() +
  theme_minimal()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=NO3)) +
  scale_color_viridis_c() +
  theme_minimal()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=NH4)) +
  scale_color_viridis_c() +
  theme_minimal()

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=IN.PO4)) +
  scale_color_viridis_c() +
  theme_minimal()


ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=Cyanobacteria)) +
  scale_color_viridis_c() +
  theme_minimal()

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




# Environmental variables can also be used with envfit which are referred to as extrinsic variables. This works best with continuous variables.If you only want to fit vector variables (continuous variables) use vectorfit and if you only want to fit factor variables (categorical variables) use factorfit but envfit can do this automatically.

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



# 3a - Pairwise beta diversity metrics!! ####
# from Jordy 

# cluster::daisy ---- Compute all the pairwise dissimilarities (distances) between observations in the data set. The original variables may be of mixed types. In that case, or whenever metric = "gower" is set, a generalization of Gower's formula is used

# “Gower's distance” is chosen by metric "gower" or automatically if some columns of x are not numeric. Also known as Gower's coefficient (1971), expressed as a dissimilarity, this implies that a particular standardisation will be applied to each variable, and the “distance” between two units is the sum of all the variable-specific distances

# The original version of daisy is fully described in chapter 1 of Kaufman and Rousseeuw (1990). Compared to dist whose input must be numeric variables, the main feature of daisy is its ability to handle other variable types as well (e.g. nominal, ordinal, (a)symmetric binary) even when different types occur in the same data set.


# calculating beta diversity for the other variables in the dataset -- basically, how different are sites based on their environmental variables? 

# create matrix for metadat
daisy_sd_data <- sd_data |>
  select(-Group, -WaterbodyName, -CollDate, -Year, -month, -Latitude, -Longitude, -Diatom, -`Green algae`, -unknown, -Crustacean, -Cyanobacteria, -Dinoflagellate) |>
  as.data.frame()
rownames(daisy_sd_data) <- sd_data$Group
daisy_sd_data  <- daisy_sd_data[,-1]
daisy_sd_data  <- as.matrix(daisy_sd_data)
daisy_sd_data  <- replace(daisy_sd_data , is.na(daisy_sd_data ), 0)
 
library(cluster)
daisy.mat <- as.matrix(daisy(scale(daisy_sd_data), metric="gower"))

daisy.mat[upper.tri(daisy.mat, diag = T)] <- NA

library(reshape2)
env_daisy<-melt(daisy.mat, na.rm=TRUE) |>
  left_join(sd_data |> select(Group,WaterbodyName, Year, month),
            by=c('Var1'='Group')) |>
  rename(Var1_WBN = WaterbodyName,
         Var1_Yr = Year,
         Var1_mon = month) |>
  left_join(sd_data |> select(Group,WaterbodyName, Year, month),
            by=c('Var2'='Group')) |>
  rename(Var2_WBN = WaterbodyName,
         Var2_Yr = Year,
         Var2_mon = month)

ggplot(env_daisy) +
  geom_histogram(aes(value))

# how different is the lake across years and months (sames sites)?
same_sites_daisy <- env_daisy |>
  filter(Var1_WBN == Var2_WBN)

ggplot(same_sites_daisy) +
  geom_point(aes(Var1_mon, value, shape=as.character(Var1_Yr), color=Var2_mon)) +
  facet_wrap(~Var1_WBN*Var2_Yr)

date_distance_daisy <- env_daisy |>
  mutate(Var1_fakedate = as.Date(paste0(Var1_Yr,'-',Var1_mon,'-01'), format='%Y-%b-%d'),
         Var2_fakedate = as.Date(paste0(Var2_Yr,'-',Var2_mon,'-01'), format='%Y-%b-%d'),
         date_distance = Var2_fakedate-Var1_fakedate)

### 3aA comparisons by temporal distance ####

ggplot(date_distance_daisy, aes(date_distance, value)) +
  geom_jitter()



### 3aB get geodistances between points ####

# get distance between sites from Jordy code: https://github.com/jvoneggers/WYLakeSedMicrobes/blob/main/VonEggers_WYMountainLakeSedimentMicrobes_DataAnalysis_Jan2024.Rmd
#Lat&lon to distance in meters
library(geosphere)
xy <- sd_data |> select(Longitude, Latitude) |> as.data.frame()
rownames(xy)<-sd_data$Group
dist_m_output<-distm(xy)/1000 # convert m to km
rownames(dist_m_output)<-rownames(xy)
names(dist_m_output)<-rownames(xy)
dist_m_output<-as.dist(dist_m_output)

#transform all distance matrices into dataframes with pairwise comparisons
coord.dist.ls<-as.matrix(dist_m_output)
coord.dist.ls[upper.tri(coord.dist.ls, diag = T)] <- NA
coord.dist.ls<-reshape2::melt(coord.dist.ls, na.rm=T)





# now look at pairwise comparisons of phyto communities using bray-curtis

# 0 is similar
comm.dist <- dist 

#pull out metadata
metadata<-sd_data[order(match(rownames(sd_data),rownames(comm.dist))),]

comm.dist.ls<-as.matrix(comm.dist)
comm.dist.ls[upper.tri(comm.dist.ls, diag = T)] <- NA
comm.dist.ls<-reshape2::melt(comm.dist.ls, na.rm=T)


combo_phyto_env_kmdist <- left_join((comm.dist.ls |> rename(phyto_dist = value)), 
                             (env_daisy |> rename(env_dist = value))) |>
  left_join(coord.dist.ls |> rename(site_dist_km = value))




### 3aC comparisons by distance!!!####
ggplot(combo_phyto_env_kmdist, aes(env_dist, phyto_dist, color=site_dist_km)) +
  geom_point() +   
  scale_color_viridis_c()

ggplot(combo_phyto_env_kmdist, aes(site_dist_km, env_dist)) +
  geom_point() 

ggplot(combo_phyto_env_kmdist, aes(site_dist_km, phyto_dist)) +
  geom_point() 




# 5. ALPHA DIVERSITY ####
#Taxonomy Line Plot
#rarefy and normalizing before calculating any alpha diversity metrics- since my samples have a variety of number of sequences or organisms

phyto_rel_abund <- phyto_data |>
  select(-WaterbodyName, -CollDate, -Latitude, -Longitude, -PO4, -NH4,-TP, -TN, -NO3, -CHLA, -Year, -month, -Diatom, -`Green algae`, -unknown, -Crustacean, -Cyanobacteria, -Dinoflagellate) |>
  pivot_longer(-Group, names_to = 'taxa', values_to = 'count') |>
  group_by(Group) |>
  mutate(rel_abund = count/sum(count, na.rm=TRUE)) |>
  ungroup() |>
  select(-count) |>
  left_join(sd_data) |>
  group_by(Group, WaterbodyName, taxa) |>
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") |>
  ungroup() |>
  group_by(WaterbodyName, taxa) #|>
#   summarise(rel_abund = mean(rel_abund, na.rm=TRUE), .groups='drop')  |>
# ungroup()  
  

taxon_pool <- phyto_rel_abund |>
  group_by(WaterbodyName, taxa) |>
  summarise(mean=mean(rel_abund,na.rm=TRUE), groups='drop') |>
  ungroup() |>
  group_by(taxa) |>
  summarize(pool = max(mean) < 3,
            mean = mean(mean),
            .groups="drop")

abundances <- inner_join(phyto_rel_abund, taxon_pool) |>
  mutate(taxa = if_else(pool, "Other", taxa)) |>
  group_by(Group, WaterbodyName, taxa) |>
  summarise(rel_abund=sum(rel_abund),
            mean = min(mean),
            .groups="drop") |>
  drop_na(taxa)
  

ggplot(abundances, aes(rel_abund, taxa, color=WaterbodyName)) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.6)) +
  theme_minimal() +
  labs(y=NULL,
       x="Relative Abundance (%)")


# 7. %cyano plot from Smith 1983 ####

# multiply N:P mass by 2.11306 to get molar
NP_cyano <- sd_data |>
  pivot_longer(cols=c(TN.TP, NO3.PO4, NH4.PO4, IN.PO4), names_to = 'NP_type')

ggplot(NP_cyano) +
  geom_point(aes(value, Cyanobacteria, color=NP_type)) +
  labs(x='N:P molar ratio', y='% Cyanobacteria abundance') +
  geom_vline(xintercept = 29) +
  theme_minimal()
  

# 8. Tributary Budget Plots ####

# monthly time interval
q_interval <- 2.628e+6 # seconds/month

pre_bud <- BoysenTribs |>
  mutate(month=month(CollDate, label=TRUE, abbr=TRUE)) |>
  filter(between(Year, 2020,2022)) |>
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN', 
                                       ShortName_Revised=='Total Ammonia as N'~'NH4', 
                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                       ShortName_Revised=='Dissolved Orthophosphate'~'PO4',
                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                       ShortName_Revised=='Discharge'~'Discharge')) |>
  drop_na(ShortName_Revised) |> # removing pH and SpC for now
  select(- ChemUnits) |> # all nutrients mg/L, Discharge L/s
  pivot_wider(names_from = ShortName_Revised, values_from = ChemValue) |>
  filter(Year < 2023) 

chems_bud <- pre_bud |>
  select(CollDate, WaterbodyName, Year, month, TN, NH4, NO3, PO4, TP) 
# get rid of rows where all are empty
chems_bud$chk_row <- rowSums(chems_bud[,5:9], na.rm=TRUE)
chems_bud <- chems_bud |>
  filter(chk_row>0) |>
  select(-chk_row)

dis_bud <- pre_bud |>
  select(CollDate, WaterbodyName, Year, month, Discharge) |>
  drop_na(Discharge)  |>
  mutate(date_minus3 = CollDate - days(3),
         date_plus3 = CollDate + days(3))

library(fuzzyjoin)

bud <- fuzzy_left_join(chems_bud, dis_bud,
                       by=c("WaterbodyName"="WaterbodyName",
                            'CollDate'="date_minus3",
                            "CollDate"="date_plus3"),
                       match_fun=list(`==`, `>`, `<`)) |>
  select(CollDate.x, WaterbodyName.y, Year.x, month.x, TN,NH4,NO3,PO4,TP,Discharge) |>
  rename(CollDate = CollDate.x,
         WaterbodyName = WaterbodyName.y,
         Year = Year.x,
         month = month.x) |>
  drop_na(WaterbodyName) |>
  # because of the fuzzy join, multiple days joined in cases where discharge was collected daily - averaging to get around that
  group_by(CollDate, WaterbodyName, Year, month) |>
  summarise(TN=mean(TN),
            NH4=mean(NH4),
            NO3=mean(NO3),
            PO4=mean(PO4),
            TP=mean(TP),
            Discharge=mean(Discharge)) |>
  ungroup() |>
  # get monthly mass loading for each parameter 
  mutate(TN_load_Mg = TN*Discharge*q_interval/1e-9,
         NH4_load_Mg = NH4*Discharge*q_interval/1e-9,
         NO3_load_Mg = NO3*Discharge*q_interval/1e-9,
         PO4_load_Mg = PO4*Discharge*q_interval/1e-9,
         TP_load_Mg = TP*Discharge*q_interval/1e-9) |>
  pivot_longer(cols=c(TN_load_Mg, 
                      NH4_load_Mg,
                      NO3_load_Mg,
                      PO4_load_Mg,
                      TP_load_Mg), values_to = 'monthly_mass_Mg', names_to = "load_type") |>
  mutate(monthly_mass_Mg = ifelse(WaterbodyName=='Wind River Outlet', monthly_mass_Mg*-1, monthly_mass_Mg))


ggplot(bud) +
  geom_bar(aes(fill = WaterbodyName, x = month, y = monthly_mass_Mg), position = "stack", stat = "identity") +
  scale_fill_viridis_d(option = "inferno", direction = -1) +
  theme_minimal() +
  labs(y = "Mass Loading (Mg)", x = "") +
  facet_wrap(~load_type*Year, scales='free_y', ncol=3) +
  #scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  scale_y_continuous(n.breaks = 5)


# 9. tributary and outlet N and P timeseries ####
ggplot(BoysenTribs |>filter(!ShortName_Revised%in%c('Discharge', 'pH', 'Conductance')), aes(CollDate,ChemValue,group=WaterbodyName, color=WaterbodyName)) +
  geom_point() +
  geom_smooth(method='gam') +
  facet_wrap(~ShortName_Revised, scales='free_y')

BoysenTribs |> 
  filter(!ShortName_Revised%in%c('Discharge', 'pH', 'Conductance')) |>
           filter(between(year(CollDate), 2020, 2022)) |>
  ggplot(aes(CollDate,ChemValue,group=WaterbodyName, color=WaterbodyName)) +
  geom_point() +
  facet_wrap(~ShortName_Revised, scales='free_y') +
  theme_minimal() +
  labs(x='', y='Concentration'~(mg~L^-1))


# 10. Boysen storage  ####
read.csv('Data/reservoir_storage_af.csv',skip=7) |>
  mutate(date = as.Date(Datetime..UTC.)) |>
  filter(between(year(date), 2020,2022)) |>
  ggplot() +
  geom_line(aes(date, Result)) +
  theme_minimal() +
  labs(x='',y='Reservoir storage (acre-ft)')


# 11. Profiles of WQ data ####
plot_profile_points <- function(param,paramname) {
  ggplot(BoysenProfile) +
    geom_point(aes(CollDate, depth_m, color=param)) +
    facet_wrap(~WaterbodyName, scales='free_y') +
    scale_color_viridis_c('') +
    scale_y_reverse() +
    theme_minimal() +
    labs(x='', y='Depth (m)', title=paramname)
}

plot_profile_points(BoysenProfile$temp_C, 'Temp')
plot_profile_points(BoysenProfile$pH, 'pH')
plot_profile_points(BoysenProfile$cond_ugL, 'SpC')
plot_profile_points(BoysenProfile$DO_mgL, 'DO')
plot_profile_points(BoysenProfile$DO_percent, 'DO %')
plot_profile_points(BoysenProfile$ORP, 'ORP')




# 12. Nutrient timeseries ####

nutrient_forms <- BoysenNutrient|>
  mutate(ShortName_Revised = case_when(ShortName_Revised=='Total Nitrogen (unfiltered)'~'TN',
                                                       ShortName_Revised=='Total Ammonia as N'~'NH4',
                                                       ShortName_Revised=='Phosphorus as P (total)'~'TP',
                                                       ShortName_Revised=='Orthophosphate as P (total)'~'PO4',
                                                       ShortName_Revised=='Nitrate plus Nitrite as N'~'NO3',
                                                       ShortName_Revised=='Chlorophyll a (phytoplankton)'~'CHLA')) |>
  mutate(form = ifelse(ShortName_Revised %in% c('TP', 'PO4'), 'P', 
                       ifelse(ShortName_Revised %in% c('TN', 'NO3', 'NH4'), 'N', NA))) |>
  mutate(BelowDet = ifelse(BelowDet==1, 'Below detection', 'Result fine'))

 
ggplot(nutrient_forms |> 
         filter(ShortName_Revised=='CHLA'),
       aes(CollDate, ChemValue, shape=BelowDet)) +
  geom_point() +
  facet_wrap(~WaterbodyName, scales='free_y') +
  scale_shape_manual('', values=c(3,16)) +
  theme_minimal() +
  labs(x='', y='Surface Chlorophyll-a Concentration'~(mu*g~L^-1))



ggplot(nutrient_forms |> 
         filter(form=='N'),
       aes(CollDate, ChemValue, shape=BelowDet, color=ShortName_Revised)) +
  geom_jitter() +
  facet_wrap(~WaterbodyName, scales='free_y') +
  scale_shape_manual('', values=c(3,16)) +
  scale_color_viridis_d('', option='plasma') +
  theme_minimal() +
  labs(x='', y='Surface Concentration'~(mg~L^-1))


ggplot(nutrient_forms |> 
         filter(form=='P'),
       aes(CollDate, ChemValue, shape=BelowDet, color=ShortName_Revised)) +
  geom_jitter() +
  facet_wrap(~WaterbodyName, scales='free_y') +
  scale_shape_manual('', values=c(3,16))  +
  scale_color_viridis_d('', option='plasma') +
  theme_minimal() +
  labs(x='', y='Surface Concentration'~(mg~L^-1))



ggplot(BoysenChem |> 
         filter(ShortName_Revised=='Secchi Depth'),
       aes(CollDate, ChemValue, color = WaterbodyName)) +
  geom_point() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='', y='Secchi depth (m)')




ggplot(BoysenChem |> 
         filter(ShortName_Revised=='Conductance'),
       aes(CollDate, ChemValue, color = WaterbodyName)) +
  geom_point() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='', y='Surface SpC'~(µS~cm^-1))



ggplot(BoysenChem |> 
         filter(ShortName_Revised=='pH'),
       aes(CollDate, ChemValue, color = WaterbodyName)) +
  geom_point() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='', y='Surface pH')


ggplot(BoysenChem |> 
         filter(ShortName_Revised=='DO, mg/L'),
       aes(CollDate, ChemValue, color = WaterbodyName)) +
  geom_point() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='', y='Surface DO'~(mg~L^-1))




# 13. storage-discharge relationships ? ####
## not good -- all have really bad R2. Data don't even look linear when plotted. This just won't work unfortunately 
storage <- read.csv('Data/reservoir_storage_af.csv',skip=7) |>
  mutate(CollDate = as.Date(Datetime..UTC.)) |>
  rename(Storage_AF = Result)

storage_discharge <- BoysenTribs |>
  filter(#WaterbodyName == 'Wind River Outlet',
         ShortName_Revised == 'Discharge') |>
  left_join(storage |> select(CollDate, Storage_AF))

ggplot(storage_discharge,aes(Storage_AF, ChemValue)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~WaterbodyName, scales='free')

WR_out_lm <- lm(ChemValue~Storage_AF, storage_discharge|>filter(WaterbodyName=='Wind River Outlet'))
summary(WR_out_lm)

WR_in_lm <- lm(ChemValue~Storage_AF, storage_discharge|>filter(WaterbodyName=='Wind River Inlet'))
summary(WR_in_lm)

MC_lm <- lm(ChemValue~Storage_AF, storage_discharge|>filter(WaterbodyName=='Muddy Creek'))
summary(MC_lm)

FMC_lm <- lm(ChemValue~Storage_AF, storage_discharge|>filter(WaterbodyName=='Fivemile Creek'))
summary(FMC_lm)




# 14. basic cyano plots ####
ggplot(sd_data, aes(month, Cyanobacteria, color=WaterbodyName, group=WaterbodyName)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='',y='% Biovolume cyanobacteria') +
  facet_wrap(~Year) 

ggplot(sd_data, aes(TN, Cyanobacteria, color=WaterbodyName, group=WaterbodyName)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='TN concentration'~(mg~L^-1),y='% Biovolume cyanobacteria') +
  facet_wrap(~Year) 

ggplot(sd_data, aes(NO3, Cyanobacteria, color=WaterbodyName, group=WaterbodyName)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='Nitrate concentration'~(mg~L^-1),y='% Biovolume cyanobacteria') +
  facet_wrap(~Year) 

ggplot(sd_data, aes(NH4, Cyanobacteria, color=WaterbodyName, group=WaterbodyName)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='Ammonium concentration'~(mg~L^-1),y='% Biovolume cyanobacteria') +
  facet_wrap(~Year) 

ggplot(sd_data, aes(TP, Cyanobacteria, color=WaterbodyName, group=WaterbodyName)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='TP concentration'~(mg~L^-1),y='% Biovolume cyanobacteria') +
  facet_wrap(~Year) 

ggplot(sd_data, aes(PO4, Cyanobacteria, color=WaterbodyName, group=WaterbodyName)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d('', option='turbo') +
  theme_minimal() +
  labs(x='Phosphate concentration'~(mg~L^-1),y='% Biovolume cyanobacteria') +
  facet_wrap(~Year) 


