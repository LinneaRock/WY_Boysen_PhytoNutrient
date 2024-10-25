#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Principle Component Analysis to look at heterogenity across reservoir sites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. Reading in data and libraries ####
source('Data/CALL_DATA_LIB.R')

library(MASS)

# 2. Format data ####
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

metadat <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, Cyanobacteria, ShortName_Revised, ChemValue) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) |>
  mutate(IN = NO3+NH4) |>
  left_join(hypo_dat)

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
# the first two explain about 44% of the variance..
EigenProp[1] + EigenProp[2] + EigenProp[3]
# the first two explain about 56% of the variance..

# not great var explained stats, but maybe useful in some way still

## evaluate loadings ####
# only keep those which have variance explained >0.3
for (i in 1:14){
  print(paste0("Principal component",i))
  print(loadings[which(abs(loadings[,i])>0.3),i])
}

# PC1: higher PO4, hypolimnion DO + lower IN:PO4, Temp, Stability
# PC2: higher TP + lower SpC, hypolimnion pH, hypolimnion SpC
# PC3: higher Secchi, stability, maxdepth + lower hypolimnion temp
# PC4: higher DO + lower IN, hypolimnion temp
# PC5: higher temp, secchi + lower TP, TN, pH


# ..... more through PC 14: 

# [1] "Principal component1"
# PO4     IN.PO4       Temp  Stability    hypo_DO 
# 0.3524271 -0.3438655 -0.4070489 -0.3421449  0.3784427 
# [1] "Principal component2"
# TP        SpC    hypo_pH   hypo_SpC 
# 0.3395137 -0.4198748 -0.3399444 -0.4005812 
# [1] "Principal component3"
# Secchi  Stability   maxdepth  hypo_temp 
# 0.3277845  0.3026838  0.5304146 -0.3626276 
# [1] "Principal component4"
# DO         IN  hypo_temp 
# 0.3277073 -0.4299712 -0.3555028 
# [1] "Principal component5"
# TP         TN         pH     Secchi          H 
# -0.3145366 -0.4471863 -0.4014573  0.3477722  0.3159931 
# [1] "Principal component6"
# PO4      TN.TP     IN.PO4          H         IN 
# 0.3085940  0.3008115 -0.3496151 -0.5511829 -0.3675495 
# [1] "Principal component7"
# TN      TN.TP   hypo_SpC 
# 0.5219078  0.5042315 -0.3458683 
# [1] "Principal component8"
# TN.TP     Secchi          H 
# 0.3567446 -0.4927619  0.4233111 
# [1] "Principal component9"
# IN.PO4        SpC    hypo_pH 
# -0.3155561 -0.4179708  0.5103013 
# [1] "Principal component10"
# H         IN 
# -0.4628964  0.6019804 
# [1] "Principal component11"
# TN.TP         DO         pH     Secchi 
# 0.3087495 -0.3864807 -0.4263721 -0.3954036 
# [1] "Principal component12"
# TP  Stability   maxdepth   hypo_SpC 
# -0.4133408 -0.3629318  0.5081991 -0.3932766 
# [1] "Principal component13"
# TP         TN      TN.TP 
# -0.6212835  0.4231918 -0.3375877 
# [1] "Principal component14"
# DO         pH       Temp  Stability    hypo_pH 
# 0.4751695 -0.4844963  0.3512732  0.3142548  0.3139925 


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
#   |(Intercept) | 67.838| 2.761| 24.573| <0.0001|
#   |PC1         | -7.035| 1.349| -5.215| <0.0001|

summary(lm(Cyanobacteria ~ PC1, test_dat))
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
#   |(Intercept) | 43.327| 0.003| 15,623.781| <0.0001|
#   |PC3         |  0.021| 0.002|     11.015| <0.0001|
#   |PC2         | -0.007| 0.001|     -5.027| <0.0001|
#   |PC1         | -0.007| 0.001|     -4.981| <0.0001|

summary(lm(Latitude~PC3, test_dat)) # r2=0.362
summary(lm(Latitude~PC2, test_dat)) # r2=0.071
summary(lm(Latitude~PC1, test_dat)) # r2=0.070

summary(lm(Latitude~PC3+PC1, test_dat))
#  r2 = 0.434 , p<2.2e-16

summary(lm(Latitude~PC3+PC2, test_dat))
#  r2 = 0.435 , p<2.2e-16

summary(lm(Latitude~PC2+PC1, test_dat))
# worst, r2 = 0.141 , p<2.2e-16




ggplot(test_dat) +
  geom_point(aes(PC1, PC3, color=Latitude))
ggplot(test_dat) +
  geom_point(aes(PC1, PC3, color=WaterbodyName))


ggplot(test_dat |>
         mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay')))) +
  geom_point(aes(PC1, PC3, fill=WaterbodyName),shape=21) +
  facet_wrap(~month) +
  theme_bw() + 
  scale_fill_viridis_d('', option = 'magma') 
ggsave('Figures/pca_fig.png', width=6.5, height=4.5, units = 'in',dpi=1200) 

EigenProp[1] + EigenProp[3]
# 1+3 explain about 35% of the variance..
