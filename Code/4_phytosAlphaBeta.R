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

# monthly discharge as proxy for nutrient loading (see 2_nutrientsReservoir&Tribs.R)
loadQ <- TribLoadFlux |>
  mutate(fakedate = paste0(Year, '-', month, '-01'),
         fakedate = as.Date(fakedate, format='%Y-%b-%d')) |>
  filter(!grepl('Outlet', WaterbodyName)) |>
  group_by(fakedate) |>
  summarise(discharge = sum(Discharge, na.rm=TRUE)) |>
  ungroup () |>
  bind_rows(BoysenTribs |>
              filter(ShortName_Revised=='Discharge',
                     Year==2023,
                     !grepl('Outlet', WaterbodyName)) |>
              mutate(month=month(CollDate)) |>
              group_by(WaterbodyName, month, Year) |>
              summarise(Discharge=mean(ChemValue)) |>
              ungroup() |>
              filter(month==9) |>
              mutate(fakedate=as.Date('2023-09-01')) |>
              group_by(fakedate) |>
              summarise(discharge=sum(Discharge)) |>
              ungroup()) |>
  mutate(Year=year(fakedate),
         month=month(fakedate, label=TRUE)) |>
  select(-fakedate)

# hypolimnion WQ data
hypo_dat <- BoysenProfile |>
  filter(depth_m==maxdepth) |>
  group_by(WaterbodyName, CollDate) |>
  summarise(hypo_temp=mean(temp_C),
            hypo_pH=mean(pH),
            hypo_SpC=mean(cond_uScm),
            hypo_DO=mean(DO_mgL)
            ) |>
  ungroup() 

# create metadata 
sd_data <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  left_join(monthly_storage) |>
  left_join(loadQ) |>
  left_join(hypo_dat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  dplyr::select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate, ave_storage_AF, discharge, hypo_temp, hypo_pH, hypo_SpC, hypo_DO) |>
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
  select(-WaterbodyName, -CollDate, -julianday, -Latitude, -Longitude, -PO4,-TP, -TN, -IN, -CHLA, -TN.TP, -IN.PO4, -Year, -month, -Diatom, -`Green algae`, -Flagellate, -`Golden algae`, -Cyanobacteria, -Dinoflagellate, -DO, -Secchi, -pH, -SpC,-H, -Stability,-Temp,-maxdepth,-ave_storage_AF,-discharge, -c(hypo_temp, hypo_pH, hypo_SpC, hypo_DO)) |>
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
  ggplot(aes(rel_abund, reorder(taxa, rel_abund), fill=WaterbodyName, color=WaterbodyName)) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.6), shape=21) +
  theme_classic() +
  labs(y=NULL,
       x="Relative Abundance (%)") +
  scale_fill_viridis_d('',option='magma') +
  scale_color_viridis_d('',option='magma') +
  theme(legend.position = c(0.8,0.35)) +
  theme(legend.position = 'none')
# ggsave('Figures/Fig4/rel_abundance.png',height=4.5,width=6.5,units='in',dpi=1200)





# 3. BETA DIVERSITY ANALYSIS ####


# create distance matrix of phytoplankton communities
dist_phyto <- phyto_data |>
  select(-WaterbodyName, -CollDate, -julianday, -Latitude, -Longitude, -PO4, -TP, -TN, -CHLA, -TN.TP,  -IN.PO4, -Year, -month, -Diatom, -`Green algae`, -Flagellate, -`Golden algae`, -Cyanobacteria, -Dinoflagellate, -DO, -Secchi, -pH, -SpC,-H, -Temp, -Stability,-maxdepth, -ave_storage_AF,-IN, -discharge, -c(hypo_temp, hypo_pH, hypo_SpC, hypo_DO))
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
adonis2(dist~Longitude, sd_data)
set.seed(69420)
adonis2(dist~WaterbodyName, sd_data) # high p-value == sites are the same in terms of their beta diversity (i.e., comparing samples to each other and answers question 'how different')? p=0.97
set.seed(69420)
adonis2(dist~month, sd_data)  # p=0.021
adonis2(dist~julianday, sd_data) 
adonis2(dist~CHLA, sd_data, na.rm=TRUE) # chla not collected 2021-05-18, permanova cannot run
adonis2(dist~TN, sd_data) 
adonis2(dist~TP, sd_data)
adonis2(dist~TN*TP, sd_data)
adonis2(dist~TN.TP, sd_data)
adonis2(dist~IN, sd_data)
set.seed(69420)
adonis2(dist~PO4, sd_data) #p=0.069
set.seed(69420)
adonis2(dist~IN.PO4, sd_data) # p=0.011
set.seed(69420)
adonis2(dist~Cyanobacteria, sd_data) #p=0.037
set.seed(69420)
adonis2(dist~Temp, sd_data) #p=0.028
adonis2(dist~hypo_temp, sd_data)
adonis2(dist~pH, sd_data)
adonis2(dist~hypo_pH, sd_data)
adonis2(dist~SpC, sd_data)
adonis2(dist~hypo_SpC, sd_data)
adonis2(dist~DO, sd_data) 
set.seed(69420)
adonis2(dist~hypo_DO, sd_data) #p=0.046
adonis2(dist~Secchi, sd_data)
adonis2(dist~H, sd_data) 
set.seed(69420)
adonis2(dist~Stability, sd_data) # sig p =0.026
adonis2(dist~maxdepth, sd_data)
set.seed(69420)
adonis2(dist~ave_storage_AF, sd_data) # p=0.044
set.seed(69420)
adonis2(dist~discharge, sd_data) # p=0.099

psych::pairs.panels(sd_data |> select(-c(1:8)))

## 3B. NMDS - beta dispersion plots ####
set.seed(69420)
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
  geom_point(aes(fill=WaterbodyName),size=2,shape=21) +
  theme_minimal() +
  scale_fill_viridis_d('', option='magma') +
  geom_text(label='dist~location \np = 0.97', mapping = aes(x = 1, y = 2)) +
  theme(legend.position = 'none')
ggsave('Figures/nmds_space.png',height=4.5,width=6.5,units='in',dpi=1200)

# however, both space and the % cyanobacteria biomass in the community were significant contributors 
# cyanos are consistent through time 
ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(fill=month),shape=21,size=2) +
  theme_minimal() +
  facet_wrap(~Year) +
  scale_fill_viridis_d('',option='plasma')

ggplot(scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(fill=Cyanobacteria),shape=21,size=2) +
  scale_fill_viridis_c() +
  theme_minimal()

## 3C. INTRINSIC VARIALBES - investigate the species which may be driving the site distribution pattern ####
set.seed(69420)
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
  scale_color_viridis_d('', option='plasma') +
  geom_segment(sig.spp.fit, mapping=aes(x=0, xend=NMDS1*2, y=0, yend=NMDS2*2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(sig.spp.fit, mapping=aes(x=NMDS1*2, y=NMDS2*2, label = spp.variables), cex = 4, direction = "both", segment.size = 0.25) + #add labels, use ggrepel::geom_text_repel so that labels do not overlap
  geom_text(label='dist~month \np = 0.021', mapping = aes(x = 1, y = 2)) +
  theme(legend.position= 'none')
# ggsave('Figures/fig4/intrinsicvariables.png',height = 4.5,width=6.5,units='in',dpi=1200)




## 3D EXTRINSIC VARIABLES ####
#-- Environmental variables can also be used with envfit which are referred to as extrinsic variables. This works best with continuous variables.If you only want to fit vector variables (continuous variables) use vectorfit and if you only want to fit factor variables (categorical variables) use factorfit but envfit can do this automatically.
set.seed(69420)
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
  geom_segment(sig.envbio.fit, mapping=aes(x=0, xend=NMDS1*2, y=0, yend=NMDS2*2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(sig.envbio.fit, mapping=aes(x=NMDS1*2, y=NMDS2*2, label = envbio.variables), cex = 4, direction = "both", segment.size = 0.25) + #add labels, use ggrepel::geom_text_repel so that labels do not overlap
  geom_text(label='dist~%cyano \np = 0.037', mapping = aes(x = 1, y = 2))  +
  theme(legend.position='none')
#ggsave('Figures/fig4/extrinsicvariables.png',height = 4.5,width=6.5,units='in',dpi=1200)


intr + ext +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')
ggsave('Figures/intrinsicextrinsicNMDS.png',height = 6.5,width=10.5,units='in',dpi=1200)

#ggsave('Figures/intrinsicextrinsic_LEGENDS.png',height = 6.5,width=10.5,units='in',dpi=1200)


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

## 4a. add toxins present to nearest sampling site ####
library(sf)
library(sp)
library(raster)
library(gdistance)
cyano_prep <- cyanotoxin |>
  mutate(month = month(CollDate,label=TRUE,abbr=TRUE),
         Year=year(CollDate)) |>
  dplyr::select(-CollDate, -Advisory) |>
  mutate(toxinpresent=as.character(toxinpresent)) |>
  mutate(toxinpresent=ifelse(toxinpresent=='N',0,1)) |>
  filter(toxinpresent != 0) |>
  distinct() |>
  mutate(ID = paste0(Year, month, row_number()))

# read in Boysen shapefile, add crs, and check units 
lake_shapefile <- st_read('C:/Users/lrock1/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')
st_crs(lake_shapefile)
lake_shapefile <- st_transform(lake_shapefile, crs = 3738)
st_crs(lake_shapefile)$units # feet, weird

# convert shapefile to raster grid
lake_raster <-raster(extent(lake_shapefile),res=100)
lake_raster <- rasterize(lake_shapefile, lake_raster, field=1)
plot(lake_raster)

# create transition layer
lake_transition <- transition(lake_raster, transitionFunction = function(x) 1/mean(x), directions = 8)
lake_transition <- geoCorrection(lake_transition, type = "c", multpl = FALSE)

# convert for working with shortest path function
cyano_locations <- cyano_prep |>
  dplyr::select(ID, Long, Lat)
coordinates(cyano_locations) <- ~Long+Lat
proj4string(cyano_locations) <- CRS("+init=epsg:4326")
cyano_locations <- spTransform(cyano_locations, CRS(st_crs(lake_shapefile)$proj4string))

# 7 sampling locations along reservoir
sample_locations_prep <- sd_data |>
  dplyr::select(Longitude, Latitude, WaterbodyName) |>
  distinct()

# convert for working with shortest path function
sample_locations <- sample_locations_prep
coordinates(sample_locations) <- ~Longitude+Latitude
proj4string(sample_locations) <- CRS("+init=epsg:4326")
sample_locations <- spTransform(sample_locations, CRS(st_crs(lake_shapefile)$proj4string))

## shortest path calculation ####
# extract coordinates
toxin_spots <- as.vector(as.data.frame(cyano_prep) |> dplyr::select(ID) |> distinct())[["ID"]]

sites <- as.vector(as.data.frame(sample_locations_prep) |> dplyr::select(WaterbodyName) |> distinct())[["WaterbodyName"]]

# create dataframe for distances
distance_df <- data.frame(ID=NA,
                          WaterbodyName=NA,
                          distance_km = NA)

# run function over all sites and build dataframe
for(t in toxin_spots) {
  for(s in sites) {
    tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
      
      # subset sites
      coord_A <- coordinates(cyano_locations[cyano_locations$ID==t, ])
      coord_B <- coordinates(sample_locations[sample_locations$WaterbodyName==s, ])
      
      
      # Calculate the shortest path distance
      shortest_path <- shortestPath(lake_transition, coord_A, coord_B, output = "SpatialLines")
      
      # # Plot the lake, points, and shortest path
      # plot(st_geometry(lake_shapefile), col = "yellow", main = "Lake with Shortest Path between Points A and B")
      # plot(shortest_path, add = TRUE, col = "blue", lwd = 2)
      # plot(cyano_locations, add = TRUE, col = "red", pch = 20)
      # plot(sample_locations, add = TRUE, col = "green", pch = 20)
      
      distance <- SpatialLinesLengths(shortest_path)*0.0003048 # convert feet to km
      
      tmp <- data.frame(ID=t,
                        WaterbodyName=s,
                        distance_km=distance)
      
      distance_df <- distance_df |>
        rbind(tmp)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
    
  }
}

# minimize distances
distance_minimized_df <- distance_df |>
  group_by(ID) |>
  mutate(minDist = min(distance_km)) |>
  ungroup() |>
  drop_na() |>
  filter(distance_km == minDist) |>
  dplyr::select(-minDist)

error_toxin_locations <- cyano_prep |>
  anti_join(distance_minimized_df) |>
  st_as_sf(coords=c('Long','Lat'), crs=4326)


# manually check locations of errors
library(ggsflabel)
ggplot() +
  geom_sf(st_geometry(lake_shapefile), mapping=aes()) +
  geom_sf(error_toxin_locations, mapping=aes()) +
  geom_sf_text_repel(error_toxin_locations, mapping=aes(label=ID))


cyano <- cyano_prep |>
  st_as_sf(coords=c('Long','Lat'), crs=4326)
ggplot() +
  geom_sf(st_geometry(lake_shapefile), mapping=aes()) +
  geom_sf(st_geometry(cyano), mapping=aes())


# create dataframe with locations to use in Random Forest
toxin_distance_mins <- cyano_prep |>
  left_join(distance_minimized_df) |>
  mutate(WaterbodyName=ifelse(ID%in%c('2023Sep25','2022Oct19'),'Lacustrine Pelagic: Dam', WaterbodyName)) |>
  mutate(WaterbodyName=ifelse(ID%in%c('2021Aug3','2021Jul2'),'Tough Creek Campground', WaterbodyName)) |>
  mutate(toxinpresent='Toxin present')


cy_ts <- ggplot(cyano_density|>
                  left_join(toxin_distance_mins) |>
                  mutate(toxinpresent=ifelse(is.na(toxinpresent), 'No toxin', toxinpresent)) |>
         mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay')))) +
  geom_point(aes(month, Cyanobacteria/1000, fill=WaterbodyName, group=WaterbodyName, shape=toxinpresent)) +
  geom_path(aes(month, Cyanobacteria/1000, color=WaterbodyName, group=WaterbodyName)) +
  scale_color_viridis_d('',option='magma') +
  scale_fill_viridis_d('',option='magma') +
  scale_shape_manual(values=c(21, 8)) +
  facet_wrap(~Year, scales='free_y') +
  theme_minimal() +
  labs(x='',y='Cyanobacteria denisty'~(~1000~cells~L^-1)) +
  theme(legend.position='none')
#ggsave('Figures/Fig4/cyano_density_raw.png', width=6.5,height=4.5,units='in',dpi=1200)

cy_ts + ab +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')
ggsave('Figures/cyanoPhyto.png', width=11.5,height=6.5,units='in',dpi=1200)

library(ggridges)
ggplot(cyano_density |>
         mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay')))) +
  geom_density_ridges(aes(Cyanobacteria/1000, WaterbodyName, fill=month), alpha = 0.5, scale = 1.5, 
                      quantile_lines = FALSE, size = 0.7, color = 'grey10') +
  scale_fill_viridis_d('', option='magma') +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme_minimal() +
  labs(y='',x='Cyanobacteria denisty'~(~1000~cells~L^-1))

