#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Principle Component Analysis to look at heterogenity across reservoir sites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. Reading in data and libraries ####
source('Data/CALL_DATA_LIB.R')

library(MASS)

# 2. Format data ####
metadat <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, Cyanobacteria, ShortName_Revised, ChemValue) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) |>
  mutate(IN = NO3+NH4)

pca_dat <- metadat |>
  select(-c(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, Cyanobacteria, NO3, NH4)) # getting rid of NO3 and NH4 to reduce redundancy 

rownames(pca_dat) <- metadat$Group

# unfortunately need to get rid of CHLA because of NAs
pca_dat <- pca_dat |>
  select(-CHLA)

# 3. Run PCA ####
PCA <- prcomp(pca_dat, scale=TRUE)
spc <- summary(PCA, loadings=T, cutoff=0)

## 3a evaluate which PCAs are best ####
# here are the eigen values, or the magnitude of direction in new feature space:
EigenValues <- (PCA[['sdev']])^2
# I squared because 'sdev' is the sqr roots of eigenvalues of covariance matrix!
EigenValues

# here are the eigen vectors, or the directionality of the new feature space:
## also we can all this the loadings (I think), or the rotated EOFs that are 
## scaled to importance of expl. variance
loadings <- PCA[['rotation']]
loadings

# This calculates the proportion of total variance explained:
EigenProp <- EigenValues / sum(EigenValues)
EigenProp

# calculate the confidence intervals
# Here, I am using N* like in the Barnes Notes
# and assuming all features are independent 
ci <- EigenProp*sqrt(2/nrow(pca_dat)) 

# put data together for plotting
screedat <- data.frame(EigenProp, ci, EOF=1:length(EigenProp))

# plot the eigenvalues as percent variance explained 
ggplot(screedat) +
  #geom_line(aes(EOF, EigenProp)) +
  geom_point(aes(EOF, EigenProp)) +
  geom_errorbar(aes(EOF, EigenProp, 
                    ymin=EigenProp-ci, ymax=EigenProp+ci)) +
  theme_classic() +
  labs(x='EOF',y='Proportion variance explained')

# check out scree plot to see which EOFs/PCs to focus on
barplot(EigenProp)
abline(h = mean(EigenProp),
       col = "red")

EigenProp[1] + EigenProp[2]
# the first two explain about 43% of the variance..
EigenProp[1] + EigenProp[2] + EigenProp[3]
# the first two explain about 56% of the variance..

# not great var explained stats, but maybe useful in some way still

## evaluate loadings ####
# only keep those which have variance explained >0.3
for (i in 1:14){
  print(paste0("Principal component",i))
  print(loadings[which(abs(loadings[,i])>0.3),i])
}

# PC1: higher PO4, TP + lower TN:TP, IN:PO4, Temp, Stability
# PC2: higher SpC, DO + lower TP, Temp
# PC3: lower Secchi, maxdepth, IN
# PC4: higher TN, pH, lower H
# PC5: lower TP, IN:PO4, DO, maxdepth + higher temp


# ..... more through PC 14: 

#   [1] "Principal component1"
# PO4         TP     IN.PO4       Temp  Stability         IN 
# 0.3803084  0.3433609 -0.3782102 -0.3015427 -0.3156025  0.3215030 
# [1] "Principal component2"
# TP        SpC         DO       Temp 
# -0.3105176  0.5300010  0.4799849 -0.3405117 
# [1] "Principal component3"
# NO3     Secchi   maxdepth         IN 
# -0.4364819 -0.4481414 -0.4688616 -0.4099500 
# [1] "Principal component4"
# TN         pH          H 
# -0.3279613 -0.6109626  0.4075194 
# [1] "Principal component5"
# NH4     IN.PO4          H 
# -0.6036331 -0.3247827 -0.4320130 
# [1] "Principal component6"
# TN        NO3  Stability   maxdepth          H 
# -0.3382248  0.3199319 -0.3369883 -0.4220135  0.3601932 
# [1] "Principal component7"
# NH4         TN     IN.PO4 
# 0.5978342 -0.3120375 -0.3009485 
# [1] "Principal component8"
# TN     TN.TP         H 
# 0.5229421 0.6799276 0.3939246 
# [1] "Principal component9"
# Secchi  Stability 
# -0.5778705  0.5045720 
# [1] "Principal component10"
# TN.TP         pH     Secchi          H 
# -0.3589684  0.4302473  0.3886359  0.5255346 
# [1] "Principal component11"
# SpC         DO     Secchi   maxdepth 
# -0.5611803  0.4659381  0.3962784 -0.3517334 
# [1] "Principal component12"
# SpC         pH  Stability   maxdepth 
# 0.5058300 -0.4616193  0.4507610 -0.3609202 
# [1] "Principal component13"
# TP         TN      TN.TP 
# 0.7496656 -0.3815696  0.3797196 
# [1] "Principal component14"
# DO       Temp  Stability   maxdepth 
# -0.4838949 -0.5879313  0.3615322 -0.4252600 
# [1] "Principal component15"
# PO4     IN.PO4 
# -0.6998212 -0.6594160 
# [1] "Principal component16"
# NO3        IN 
# -0.657133  0.706793 


# reformat the data for better plotting
loadings.p <- as.data.frame(loadings)
loadings.p <- tibble::rownames_to_column(loadings.p, 'var') 
loadings.p <- loadings.p |>
  pivot_longer(c(2:15), names_to = 'EOF', values_to = 'Loading') |>
  mutate(EOF = gsub('PC', 'EOF', EOF))

ggplot(loadings.p |> filter(EOF %in% c('EOF1', 'EOF2', 'EOF3'))) +
  geom_bar(aes(var, Loading), stat='identity') +
  facet_grid(~EOF) +
  theme_classic() +
  labs(x='') +
  theme(axis.text.x = element_text(angle=45, hjust=1))


### plot the loadings for the first 3 EOFs ####
# regressions for the first 2 PCs
# pull out the scores
PC1 <- PCA[['x']][,1]
PC2 <- PCA[['x']][,2]
PC3 <- PCA[['x']][,3]
# reformat for better plotting
# scale the data for plotting because otherwise it is impossible to see any 
# patterns
pca_dat.reg <- scale(pca_dat) |>
  bind_cols(PC1=PC1)|>
  bind_cols(PC2=PC2) |>
  bind_cols(PC3=PC3) |>
  # choose the vars with >3 var explained from PC1 & PC2
  select(PC1, PC2, PC3, PO4, TP, IN, IN.PO4, Temp, Stability, SpC, DO, Secchi, maxdepth) |>
  pivot_longer(c(PO4, TP, IN, IN.PO4, Temp, Stability, SpC, DO,  Secchi, maxdepth), names_to = 'var', values_to = 'Scaled data') |>
  pivot_longer(c('PC1','PC2','PC3'), names_to = 'PC', values_to = 'PC values') |>
  mutate(PC=case_when(PC=='PC1'~paste0('PC1 (', round(EigenProp[1] * 100,1),
                                     '% variance explained)'),
                      PC=='PC2'~paste0('PC2 (', round(EigenProp[2] * 100 ,1),
                          '% variance explained)'),
                      PC=='PC3'~paste0('PC3 (', round(EigenProp[3] * 100 ,1),
                          '% variance explained)')))

ggplot(pca_dat.reg) +
  geom_point(aes(`PC values`, `Scaled data`, color=var, group=PC)) +
  facet_wrap(~PC) +
  labs(x = '') +
  scale_color_viridis_d('') +
  theme_classic()

# tranform original data into the new PC feature space and just keep PCs 1 and 2
newdataA <- as.data.frame(predict(PCA, pca_dat)[,1:2])

library(ggbiplot)
ggbiplot(PCA, labels=rownames(pca_dat), obs.scale = 1, var.scale = 1,
         labels.size = 0.5, varname.size = 5, varname.adjust = 3) +
  theme_classic() +
  geom_point(newdataA, mapping=aes(PC1,PC2), alpha=0.25)

# 4. Plots of PCs with variables of interest ####
# format data
pc.scores <- as.data.frame(PCA[["x"]])

test_dat <- metadat |>
  bind_cols(pc.scores |> select(PC1, PC2, PC3))

ggplot(test_dat) +
  geom_boxplot(aes(WaterbodyName, PC1))

ggplot(test_dat) +
  geom_boxplot(aes(month, PC1))

# check linear models - 
library(knitr)
library(broom)


# cyanos 
Cyano.model<-lm(Cyanobacteria ~ 1, data=test_dat)
# perform forward selection using BIC (with k=log(n))
step.modelBIC<-MASS::stepAIC(Cyano.model, direction="forward",
                             scope= ~ PC1 + PC2 + PC3,
                             trace=F, k=log(80931))

summary(step.modelBIC)
# display the summary table for the model chosen by forward BIC selection
tidy(step.modelBIC) %>%
  mutate(p.value = scales::pvalue(p.value, accuracy = 0.0001)) %>% kable(
    caption = "Estimates for the MLR model selected 
                through forward selection using BIC as criterion.",
    col.names = c("Predictor", "B", "SE", "t", "p"),
    digits = c(3, 3, 3, 3, 4),
    align = c("l", "r", "r", "r", "r"),
    format.args = list(big.mark = ",")
  )
# Table: Estimates for the MLR model selected 
# through forward selection using BIC as criterion.
# 
# |Predictor   |      B|    SE|      t|       p|
#   |:-----------|------:|-----:|------:|-------:|
#   |(Intercept) | 67.838| 2.872| 23.618| <0.0001|
#   |PC2         | -6.688| 1.871| -3.575|  0.0005|

summary(lm(Cyanobacteria ~ PC2, test_dat))
# pretty bad model


# space

lat.model<-lm(Latitude ~ 1, data=test_dat)

# perform forward selection using BIC (with k=log(n))
step.modelBIC<-MASS::stepAIC(lat.model, direction="forward",
                             scope= ~ PC1 + PC2 + PC3,
                             trace=F, k=log(80931))

summary(step.modelBIC)
# display the summary table for the model chosen by forward BIC selection
tidy(step.modelBIC) %>%
  mutate(p.value = scales::pvalue(p.value, accuracy = 0.0001)) %>% kable(
    caption = "Estimates for the MLR model selected 
                through forward selection using BIC as criterion.",
    col.names = c("Predictor", "B", "SE", "t", "p"),
    digits = c(3, 3, 3, 3, 4),
    align = c("l", "r", "r", "r", "r"),
    format.args = list(big.mark = ",")
  )
# Table: Estimates for the MLR model selected 
# through forward selection using BIC as criterion.
# 
# |Predictor   |      B|    SE|          t|       p|
#   |:-----------|------:|-----:|----------:|-------:|
#   |(Intercept) | 43.327| 0.003| 15,007.953| <0.0001|
#   |PC3         | -0.022| 0.002|    -11.194| <0.0001|
#   |PC1         | -0.006| 0.001|     -4.414| <0.0001|

summary(lm(Latitude~PC3+PC1, test_dat))
# not great - but could be worse, r2 = 0.5169, p<2.2e-16

ggplot(test_dat) +
  geom_point(aes(PC1, PC3, color=Latitude))
ggplot(test_dat) +
  geom_point(aes(PC1, PC3, color=WaterbodyName))


ggplot(test_dat |>
         mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay')))) +
  geom_point(aes(PC1, PC3, color=WaterbodyName)) +
  facet_wrap(~month) +
  theme_bw() + 
  scale_color_viridis_d('', option = 'turbo')
ggsave('Figures/PCA/pca_fig.png', width=6.5, height=4.5, units = 'in',dpi=1200) 

EigenProp[1] + EigenProp[3]
# the first two explain about 39% of the variance..
