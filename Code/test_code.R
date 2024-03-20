#-----------------------------------------------#
# copying Kelsey's code
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. Prep data for plotting ####
head(BoysenChemPhys)
head(BoysenPhyto) 
head(BoysenTribs)

phyto_class <- read.csv('Data/phyto_class.csv') |>
  rename(cat=Class) 

phyto_class$cat <- sub('^$', 'unknown', phyto_class$cat)



BoysenPhyto_A <- BoysenPhyto |>
  left_join(phyto_class)  |>
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



BoysenPhyto_cat <- BoysenPhyto |>
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





ggplot(BoysenPhyto_A) +
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
phyto_data <- BoysenPhyto_A |>
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



# 3. NMDS and 4. BETA DIVERSITY DISPERSION ####
#beta Dispersion plot --- I think these are the NMDS plots?????? 
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
            .groups="drop")
  

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


# 10.   
































