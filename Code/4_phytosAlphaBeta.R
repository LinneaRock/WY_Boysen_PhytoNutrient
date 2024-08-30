#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Figures 4-7, Phytoplankton community dynamics 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')
library(vegan)

# monthly reservoir storage
monthly_storage <- read.csv('Data/reservoir_storage_af.csv',skip=7) |>
  mutate(date = as.Date(Datetime..UTC.)) |>
  filter(between(year(date), 2020,2023)) |>
  mutate(month = month(date, label=TRUE),
         Year= year(date)) |>
  group_by(month, Year) |>
  summarise(ave_storage_AF = mean(Result))

# create metadata 
sd_data <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  left_join(monthly_storage) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  dplyr::select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate, ave_storage_AF) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) |>
  mutate(IN=NO3+NH4) |>
  select(-c(NO3, NH4))# getting rid of NO3 and NH4 to reduce redundancy 


# match up data for later
phyto_data <- BoysenPhyto_A |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, Genus.Species.Variety, indsum, WaterbodyName, CollDate, Year, month) |>
  left_join(sd_data) |>
  pivot_wider(names_from = Genus.Species.Variety, values_from = indsum) |>
  as.data.frame() 

# 2. ALPHA DIVERSITY ANALYSIS ####
phyto_rel_abund <- phyto_data |>
  select(-WaterbodyName, -CollDate, -julianday, -Latitude, -Longitude, -PO4,-TP, -TN, -IN, -CHLA, -TN.TP, -IN.PO4, -Year, -month, -Diatom, -`Green algae`, -Flagellate, -`Golden algae`, -Cyanobacteria, -Dinoflagellate, -DO, -Secchi, -pH, -SpC,-H, -Stability,-Temp,-maxdepth,-ave_storage_AF) |>
  pivot_longer(-Group, names_to = 'taxa', values_to = 'count') |>
  group_by(Group) |>
  mutate(rel_abund = count/sum(count, na.rm=TRUE)) |>
  ungroup() |>
  select(-count) |>
  left_join(sd_data) |>
  group_by(Group, WaterbodyName, taxa) |>
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") |>
  ungroup() |>
  #mutate to make prettier yaxis in plot below 
  mutate(taxa=ifelse(taxa=='Undetermined Pennate','Pennate', taxa)) 


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

#Taxonomy Line Plot - aggregated over all time within each samping location
ab<- abundances |>
  drop_na()|>
  arrange(desc(mean)) |>
  mutate(taxa=factor(taxa)) |>
  filter(mean>5) |> # annoying filter just to make plot less complicated by removing the lowest taxa
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  ggplot(aes(rel_abund, reorder(taxa, rel_abund), color=WaterbodyName)) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.6)) +
  theme_classic() +
  labs(y=NULL,
       x="Relative Abundance (%)") +
  scale_color_viridis_d('',option='turbo') +
  theme(legend.position = c(0.8,0.35))
# ggsave('Figures/Fig4/rel_abundance.png',height=4.5,width=6.5,units='in',dpi=1200)





# 3. BETA DIVERSITY ANALYSIS ####


# create distance matrix of phytoplankton communities
dist_phyto <- phyto_data |>
  select(-WaterbodyName, -CollDate, -julianday, -Latitude, -Longitude, -PO4, -TP, -TN, -CHLA, -TN.TP,  -IN.PO4, -Year, -month, -Diatom, -`Green algae`, -Flagellate, -`Golden algae`, -Cyanobacteria, -Dinoflagellate, -DO, -Secchi, -pH, -SpC,-H, -Temp, -Stability,-maxdepth, -ave_storage_AF,-IN)
rownames(dist_phyto) <- dist_phyto$Group
dist_phyto <- dist_phyto[,-1]
dist_phyto <- as.matrix(dist_phyto)
dist_phyto <- replace(dist_phyto, is.na(dist_phyto), 0)


dist <- vegdist(dist_phyto, method = 'bray')

dist

## 3A. PERMANOVA with adonis2 can help test relationships between water quality values (for example) and community composition ####

# PERMANOVA Explained by Chat GPT
# What is PERMANOVA?
#   
#   PERMANOVA is a statistical test used to compare the differences between groups of multivariate data.
# It is particularly useful when you have complex data with many variables and you want to see if the composition of these variables differs significantly between groups.
# How Does It Work?
#   
#   Distance Matrix: First, it calculates a distance matrix that quantifies the dissimilarity between all pairs of samples based on their multivariate data. In your case, this could be a measure of how different the phytoplankton communities are from each other.
# Permutations: It then permutes (randomly rearranges) the data many times to create a distribution of possible outcomes under the null hypothesis (that there are no differences between groups).
# Comparison: It compares the observed differences between groups to this distribution to determine if the observed differences are statistically significant.
# Why Use PERMANOVA?
#   
#   Non-parametric: Unlike traditional ANOVA, PERMANOVA does not assume normal distribution of the data, making it suitable for ecological data, which often do not meet these assumptions.
# Multivariate: It can handle multiple variables at once, making it ideal for complex datasets like those involving multiple species of phytoplankton.

adonis2(dist~Latitude, sd_data)
set.seed(06261993)
adonis2(dist~WaterbodyName, sd_data) # high p-value == sites are the same in terms of their beta diversity (i.e., comparing samples to each other and answers question 'how different')? p=0.97
set.seed(06261993)
adonis2(dist~month, sd_data)  # p=0.024
adonis2(dist~julianday, sd_data) 
adonis2(dist~CHLA, sd_data, na.rm=TRUE) # chla not collected 2021-05-18, permanova cannot run
adonis2(dist~TN, sd_data) 
adonis2(dist~TP, sd_data)
adonis2(dist~TN*TP, sd_data)
adonis2(dist~TN.TP, sd_data)
adonis2(dist~IN, sd_data)
set.seed(06261993)
adonis2(dist~PO4, sd_data) #p=0.098
set.seed(06261993)
adonis2(dist~IN.PO4, sd_data) # p=0.02
set.seed(06261993)
adonis2(dist~Cyanobacteria, sd_data) #p=0.036
adonis2(dist~pH, sd_data)
adonis2(dist~SpC, sd_data)
adonis2(dist~DO, sd_data) 
adonis2(dist~Secchi, sd_data)
adonis2(dist~H, sd_data) 
set.seed(06261993)
adonis2(dist~Stability, sd_data) # sig p =0.017
adonis2(dist~maxdepth, sd_data)
set.seed(06261993)
adonis2(dist~ave_storage_AF, sd_data) # p=0.058

psych::pairs.panels(sd_data |> select(-c(1:8)))

## 3B. NMDS - beta dispersion plots ####
set.seed(06261993)
nmds <- metaMDS(dist)
nmds # make note of the stress value, this shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good -- 0.18

# pull out scores for beautiful plotting
scores <- scores(nmds) |>
  as_tibble(rownames='Group') |>
  # and join to metadata
  left_join(sd_data)

# plot showing surprising result of no effect of space on phytoplankton communities
scores |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=WaterbodyName),size=2) +
  theme_minimal() +
  scale_color_viridis_d('', option='turbo') +
  geom_text(label='dist~location \np = 0.97', mapping = aes(x = 1, y = 2)) 
ggsave('Figures/nmds_space.png',height=4.5,width=6.5,units='in',dpi=1200)

# however, both space and the % cyanobacteria biomass in the community were significant contributors 
# cyanos are consistent through time 
ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(fill=month),shape=21,size=2) +
  theme_minimal() +
  facet_wrap(~Year) +
  scale_fill_viridis_d('',option='magma')

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(fill=Cyanobacteria),shape=21,size=2) +
  scale_fill_viridis_c() +
  theme_minimal()

## 3C. INTRINSIC VARIALBES - investigate the species which may be driving the site distribution pattern ####
set.seed(06261993)
spp.fit <- envfit(nmds, dist_phyto, permutations=999)
head(spp.fit)

spp.fit_df <- as.data.frame(scores(spp.fit, display='vectors'))  #extracts relevant scores from sppifit
spp.fit_df <- cbind(spp.fit_df, spp.variables = rownames(spp.fit_df)) #and then gives them their names

spp.fit_df <- cbind(spp.fit_df, pval = spp.fit$vectors$pvals) # add pvalues to dataframe
sig.spp.fit <- subset(spp.fit_df, pval<=0.001) #subset data to show variables significant at 0.001 - there are 24 phytos that significantly explain the beta diversity of the communities below alpha=0.05. So we narrowed down to the top 8 for visualisation, which are significant below alpha=0.001

sig.spp.fit

#Now we have the relevant information for plotting the ordination!
intr<-ggplot() +
  geom_point(scores, mapping=aes(x=NMDS1, y=NMDS2, color=month)) +
  theme_minimal() +
  scale_color_viridis_d('', option='magma') +
  geom_segment(sig.spp.fit, mapping=aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(sig.spp.fit, mapping=aes(x=NMDS1, y=NMDS2, label = spp.variables), cex = 3, direction = "both", segment.size = 0.25) + #add labels, use ggrepel::geom_text_repel so that labels do not overlap
  geom_text(label='dist~month \np = 0.024', mapping = aes(x = 1, y = 2))
# ggsave('Figures/fig4/intrinsicvariables.png',height = 4.5,width=6.5,units='in',dpi=1200)




## 3D EXTRINSIC VARIABLES ####
#-- Environmental variables can also be used with envfit which are referred to as extrinsic variables. This works best with continuous variables.If you only want to fit vector variables (continuous variables) use vectorfit and if you only want to fit factor variables (categorical variables) use factorfit but envfit can do this automatically.
set.seed(06261993)
envbio.fit <- envfit(nmds, sd_data |> # get rid of these phyto groups here
                       select(-Diatom, -`Green algae`, -Flagellate, -`Golden algae`, -Cyanobacteria, -Dinoflagellate, - julianday) |>
                       # some renaming for nicer looking plots
                       rename(Storage = ave_storage_AF,
                              `IN:PO4`=IN.PO4), permutations=999, na.rm=TRUE)
head(envbio.fit)

envbio.fit_df <- as.data.frame(scores(envbio.fit, display='vectors'))  #extracts relevant scores from envifit
envbio.fit_df <- cbind(envbio.fit_df, envbio.variables = rownames(envbio.fit_df)) #and then gives them their names

envbio.fit_df <- cbind(envbio.fit_df, pval = envbio.fit$vectors$pvals) # add pvalues to dataframe
sig.envbio.fit <- subset(envbio.fit_df, pval<=0.05) #subset data to show variables significant at 0.05

sig.envbio.fit

#Now we have the relevant information for plotting the ordination!

ext<-ggplot() +
  geom_point(scores, mapping=aes(x=NMDS1, y=NMDS2, color=Cyanobacteria)) +
  theme_minimal() +
  scale_color_viridis_c('% Cyanobacteria of \ncommunity biomass') +
  theme_minimal() +
  geom_segment(sig.envbio.fit, mapping=aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(sig.envbio.fit, mapping=aes(x=NMDS1, y=NMDS2, label = envbio.variables), cex = 3, direction = "both", segment.size = 0.25) + #add labels, use ggrepel::geom_text_repel so that labels do not overlap
  geom_text(label='dist~%cyano \np = 0.036', mapping = aes(x = 1, y = 2)) 
#ggsave('Figures/fig4/extrinsicvariables.png',height = 4.5,width=6.5,units='in',dpi=1200)


intr + ext +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')
ggsave('Figures/intrinsicextrinsicNMDS.png',height = 6.5,width=10.5,units='in',dpi=1200)


# 4. Cyano biomass timeseries ####
cyano_density <- BoysenPhyto |>
  left_join(phyto_class) |>
  filter(RepNum == 0) |>
  distinct() |>
  group_by(WaterbodyName, CollDate, month, Year, cat) |>
  summarise(Density_cellsL=sum(`Density (cells/L)`)) |>
  ungroup() |>
  distinct() |>
  pivot_wider(names_from = cat, values_from = Density_cellsL) 

cyano_density <- replace(cyano_density, is.na(cyano_density), 0)

cyano_density_summarise <- cyano_density |>
  group_by(CollDate,month, Year, WaterbodyName, Cyanobacteria) |>
  summarise(`Other Phytoplankton` = sum(Diatom, `Green algae`, Dinoflagellate, `Golden algae`, Flagellate)) |>
  pivot_longer(c(Cyanobacteria, `Other Phytoplankton`), names_to = 'type', values_to = 'density') |>
  group_by(CollDate,month, Year, type) |>
  summarise(mean=mean(density),
            min=min(density),
            max=max(density)) |>
  ungroup()



ggplot() +
  geom_point(cyano_density_summarise, mapping=aes(month, mean/1000, width=0.2, color=type),
             position = position_dodge(width = 0.75)) +
  geom_errorbar(cyano_density_summarise,mapping=aes(month, mean/1000, ymin=min/1000, ymax=max/1000, width=0.2, color=type),
                position = position_dodge(width = 0.75)) +
  theme_minimal() +
  facet_wrap(~Year, scales='free_y') +
  labs(x='',y='Density'~(~1000~cells~L^-1))+
  scale_color_manual('',values=c('#999933', '#44AA99'))

ggsave('Figures/cyano_density_mean.png', width=6.5,height=4.5,units='in',dpi=1200)

cy_ts <- ggplot(cyano_density|>
         mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay')))) +
  geom_point(aes(month, Cyanobacteria/1000, color=WaterbodyName, group=WaterbodyName)) +
  geom_path(aes(month, Cyanobacteria/1000, color=WaterbodyName, group=WaterbodyName)) +
  scale_color_viridis_d('',option='turbo') +
  facet_wrap(~Year, scales='free_y') +
  theme_minimal() +
  labs(x='',y='Cyanobacteria denisty'~(~1000~cells~L^-1)) +
  theme(legend.position='none')
#ggsave('Figures/Fig4/cyano_density_raw.png', width=6.5,height=4.5,units='in',dpi=1200)

cy_ts + ab +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')
ggsave('Figures/cyanoPhyto.png', width=11.5,height=4.5,units='in',dpi=1200)

library(ggridges)
ggplot(cyano_density |>
         mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay')))) +
  geom_density_ridges(aes(Cyanobacteria/1000, WaterbodyName, fill=month), alpha = 0.5, scale = 1.5, 
                      quantile_lines = FALSE, size = 0.7, color = 'grey10') +
  scale_fill_viridis_d('', option='magma') +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme_minimal() +
  labs(y='',x='Cyanobacteria denisty'~(~1000~cells~L^-1))

