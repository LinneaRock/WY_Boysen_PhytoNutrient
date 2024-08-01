#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Using Random Forest to understand what can predict before/during/after toxin production in Boysen Reservoir
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. Reading in data and libraries ####
source('Data/CALL_DATA_LIB.R')
library(sf)
library(sp)
library(raster)
library(gdistance)
library(randomForest)
library(rsample)
library(caret)
library(pROC)

# format subset of data 
sd_data <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) |>
  mutate(IN = NO3+NH4)

# 2. Find nearest sampling locations to locations where cyanotoxins are confirmed ####
## format Boysen shapefile ####
# read in Boysen shapefile, add crs, and check units 
lake_shapefile <- st_read('C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')
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

## define points ####
### cyanotoxin locations ####
cyano_prep <- cyanotoxin |>
  mutate(month = month(CollDate,label=TRUE,abbr=TRUE),
         Year=year(CollDate)) |>
  dplyr::select(-CollDate, -Advisory) |>
  mutate(toxinpresent=as.character(toxinpresent)) |>
  mutate(toxinpresent=ifelse(toxinpresent=='N',0,1)) |>
  filter(toxinpresent != 0) |>
  distinct() |>
  mutate(ID = paste0(Year, month, row_number()))

# convert for working with shortest path function
cyano_locations <- cyano_prep |>
  dplyr::select(ID, Long, Lat)
coordinates(cyano_locations) <- ~Long+Lat
proj4string(cyano_locations) <- CRS("+init=epsg:4326")
cyano_locations <- spTransform(cyano_locations, CRS(st_crs(lake_shapefile)$proj4string))

### 7 sampling locations along reservoir ####
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
  mutate(WaterbodyName=ifelse(ID%in%c('2021Aug3','2021Jul2'),'Tough Creek Campground', WaterbodyName))

# 3. Random Forest - categorical ####
# Notes:
# RF can't handle NAs, so let's remove the CHLA. Surprisingly, it hasn't been a significant predictor of anything yet anyway in previous model explorations

# data do not need to be normal for RF - and scaling/transforming data can be harmful

## create the dataframe for RF ####
rf.data <- sd_data  |>
  # we only want to keep data for RF that corresponds to the location and year a bloom occurred
  right_join(toxin_distance_mins |> dplyr::select(WaterbodyName, Year)) |>
  left_join(toxin_distance_mins) |>
  distinct() |>
  mutate(toxinpresent=ifelse(toxinpresent==1,'During',NA)) |>
  # i wish i could figure out how to code this instead of manually, but idk how
  mutate(toxinpresent=ifelse(WaterbodyName=='Cottonwood Creek Bay' &
                               is.na(toxinpresent),'Before',toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Fremont Bay' &
                               is.na(toxinpresent) &
                               month %in% c('May','Jun'),'Before', toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Fremont Bay' &
                               is.na(toxinpresent), 'After',toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Lacustrine Pelagic: Dam' &
                               is.na(toxinpresent) &
                               month %in% c('May','Jun', 'Jul'), 'Before', toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Lacustrine Pelagic: Dam' &
                               is.na(toxinpresent), 'After', toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Tough Creek Campground' &
                               is.na(toxinpresent), 'Before', toxinpresent)) |>
  # get rid of variables we don't need
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, -Cyanobacteria, -CHLA,-WaterbodyName,-Latitude,-Longitude, - Lat, - Long, -ID, -distance_km) |>
  distinct() |>
  mutate(toxinpresent=as.factor(toxinpresent))

# to predict before/during/after, we are including the following variables: PO4, NH4, TP,TN, NO3, IN, TN:TP, IN:PO4, SpC, DO, pH, Temp, Secchi depth, stability, maxdepth at sitethat date, shannon-wiener index

## set up RF and find best parameters ####
set.seed(06261993)
split <- initial_split(rf.data, prop=0.80, strat=toxinpresent) # add strata argument to ensure random sampling is representative of the imbalanced data
training.dat <- training(split) |> mutate_if(is.numeric, round, digits=2)
testing.dat <- testing(split) |> mutate_if(is.numeric, round, digits=2) 

# model with default parameters
rf_default <- train(toxinpresent ~.,
                    training.dat,
                    #metric='RMSE',
                    method='rf',
                    tuneGrid=expand.grid(.mtry=ncol(training.dat)/3),
                    ntree=500,
                    trControl=trainControl(method='cv', number=10))

rf_default #accuracy=   0.655
# Random Forest 
# 
# 42 samples
# 16 predictors
# 3 classes: 'After', 'Before', 'During' 
# 
# No pre-processing
# Resampling: Cross-Validated (10 fold) 
# Summary of sample sizes: 38, 37, 37, 38, 39, 37, ... 
# Resampling results:
#   
#   Accuracy  Kappa    
# 0.655     0.3772727
# 
# Tuning parameter 'mtry' was held constant at a value of 5.666667

# find best mtry
set.seed(06261993)
rf_mtry <- train(toxinpresent~.,
                 data = training.dat,
                 method = "rf",
                 # metric = "RMSE",
                 tuneGrid = expand.grid(.mtry = c(1: 10)),
                 trControl = trainControl(method = "cv",
                                          number = 10,
                                          search = "grid"))

print(rf_mtry) 
plot(rf_mtry) 
# mtry=8

# Random Forest 
# 
# 42 samples
# 16 predictors
# 3 classes: 'After', 'Before', 'During' 
# 
# No pre-processing
# Resampling: Cross-Validated (10 fold) 
# Summary of sample sizes: 37, 37, 37, 38, 38, 38, ... 
# Resampling results across tuning parameters:
#   
#   mtry  Accuracy  Kappa    
# 1    0.590     0.2433333
# 2    0.610     0.2766667
# 3    0.590     0.2433333
# 4    0.615     0.2933333
# 5    0.635     0.3266667
# 6    0.635     0.3266667
# 7    0.610     0.2766667
# 8    0.660     0.3722222
# 9    0.660     0.3766667
# 10    0.610     0.2766667
# 
# Accuracy was used to select the optimal model using the largest value.
# The final value used for the model was mtry = 8.


# find best ntrees
store_maxtrees <- list()
for (ntree in c(100, 150, 250, 300, 350, 400, 450, 500, 800, 1000, 2000)) {
  set.seed(06261993)
  rf_maxtrees <- train(toxinpresent~.,
                       data = training.dat,
                       method = "rf",
                       #metric = "RMSE",
                       tuneGrid = expand.grid(.mtry = c(1: 10)),
                       trControl = trainControl(method = "cv",
                                                number = 10,
                                                search = "grid"),
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)

summary(results_tree) # 150

# Call:
#   summary.resamples(object = results_tree)
# 
# Models: 100, 150, 250, 300, 350, 400, 450, 500, 800, 1000, 2000 
# Number of resamples: 10 
# 
# Accuracy 
# Min. 1st Qu. Median  Mean 3rd Qu. Max. NA's
# 100   0.4   0.525  0.675 0.660    0.75    1    0
# 150   0.4   0.500  0.675 0.665    0.75    1    0
# 250   0.4   0.525  0.675 0.660    0.75    1    0
# 300   0.4   0.525  0.675 0.660    0.75    1    0
# 350   0.4   0.525  0.675 0.660    0.75    1    0
# 400   0.4   0.525  0.675 0.660    0.75    1    0
# 450   0.4   0.525  0.675 0.660    0.75    1    0
# 500   0.4   0.525  0.675 0.660    0.75    1    0
# 800   0.4   0.525  0.675 0.660    0.75    1    0
# 1000  0.4   0.500  0.600 0.635    0.75    1    0
# 2000  0.4   0.500  0.600 0.635    0.75    1    0
# 
# Kappa 
#      Min.    1st Qu.    Median      Mean   3rd Qu. Max. NA's
# 100     0 0.08333333 0.4166667 0.3766667 0.5000000    1    0
# 150     0 0.00000000 0.4166667 0.3888889 0.5416667    1    0
# 250     0 0.08333333 0.4166667 0.3766667 0.5000000    1    0
# 300     0 0.08333333 0.4166667 0.3766667 0.5000000    1    0
# 350     0 0.08333333 0.4166667 0.3722222 0.5000000    1    0
# 400     0 0.08333333 0.4166667 0.3722222 0.5000000    1    0
# 450     0 0.08333333 0.4166667 0.3766667 0.5000000    1    0
# 500     0 0.08333333 0.4166667 0.3722222 0.5000000    1    0
# 800     0 0.08333333 0.4166667 0.3722222 0.5000000    1    0
# 1000    0 0.00000000 0.3333333 0.3266667 0.5000000    1    0
# 2000    0 0.00000000 0.3333333 0.3266667 0.5000000    1    0

## run best model ####
set.seed(06261993)
# fit the model with the best hyperparameters
bda_fit_rf <- randomForest(toxinpresent~.,
                           training.dat,
                           method = "rf",
                           #metric = "RMSE",
                           tuneGrid = expand.grid(.mtry = c(1: 10)),
                           trControl = trainControl(method = "cv",
                                                    number = 10),
                           importance = TRUE,
                           mtry = 8,
                           ntree = 150)

# get predicted values
testing.dat$prediction <- predict(bda_fit_rf, testing.dat)

bda_fit_rf
# Call:
#   randomForest(formula = toxinpresent ~ ., data = training.dat,      method = "rf", tuneGrid = expand.grid(.mtry = c(1:10)), trControl = trainControl(method = "cv",          number = 10), importance = TRUE, mtry = 8, ntree = 150) 
# Type of random forest: classification
# Number of trees: 150
# No. of variables tried at each split: 8
# 
# OOB estimate of  error rate: 42.86%
# Confusion matrix:
#   After Before During class.error
# After      0      0      4   1.0000000
# Before     0     14      6   0.3000000
# During     0      8     10   0.4444444

# 4. Variable importance plots ####
varImpPlot(bda_fit_rf)
plot(bda_fit_rf)

varImpPlot(bda_fit_rf, type = 1, scale = TRUE,
           n.var = ncol(rf.data) - 1, cex = 0.8,
           main = "Variable importance")

bda_fit_rf$importance


reprtree::plot.getTree(bda_fit_rf)


plt.dat<-varImpPlot(bda_fit_rf, type = 1, scale = TRUE,
                    n.var = ncol(rf.data) - 1, cex = 0.8,
                    main = "Variable importance") |>
  as.data.frame() |>
  arrange(desc(MeanDecreaseAccuracy))

plt.dat<-cbind(rownames(plt.dat), plt.dat) |>
  rename(var=1)

ggplot(plt.dat) +
  geom_bar(aes(x=MeanDecreaseAccuracy,y=reorder(var, MeanDecreaseAccuracy)), stat='identity') +
  theme_classic() +
  labs(y='',x='Mean decrease accuracy', title='Variable Importance')
ggsave('Figures/RandomForest/Var_importance_bda.png', height=4.5, width=6.5, dpi=1200)


# 5. ROC and confustion matrix ####
## plot ROC curve ####
# predict test set, get probs instead of response
predictions <- as.data.frame(predict(bda_fit_rf, testing.dat, type = "prob"))

# predict class and then attach test class
predictions$predict <- names(predictions)[1:3][apply(predictions[,1:3], 1, which.max)]
predictions$observed <- testing.dat$toxinpresent
head(predictions)


roc.Before <- roc(ifelse(predictions$observed=="Before", "Before", "non-Before"), as.numeric(predictions$Before))


# others
roc.During <- roc(ifelse(predictions$observed=="During", "During", "non-During"), as.numeric(predictions$Before))
roc.After <- roc(ifelse(predictions$observed=="After", "After", "non-After"), as.numeric(predictions$Before))



roc.dat <- data.frame(toxinpresent='Before',
                      Sensitivity=roc.Before[["sensitivities"]],
                      Specificity=roc.Before[["specificities"]]) |>
  rbind(data.frame(toxinpresent='During',
                   Sensitivity=roc.During[["sensitivities"]],
                   Specificity=roc.During[["specificities"]])) |>
  rbind(data.frame(toxinpresent='After',
                   Sensitivity=roc.After[["sensitivities"]],
                   Specificity=roc.After[["specificities"]])) |>
  mutate(toxinpresent=factor(toxinpresent, levels=c('Before','During','After')))

AUC <- data.frame(toxinpresent=c('Before','During','After'),
                  auc=c(roc.Before[["auc"]],roc.During[["auc"]],roc.After[["auc"]]))



rocplot<-ggplot() +
  theme_classic() +
  geom_path(roc.dat, mapping=aes(x=1-Specificity,y=Sensitivity,color = toxinpresent)) +
  scale_color_manual('Cyanotoxin presence', values=c('#88CCEE','#999933','#44AA99')) +
  geom_abline(slope=1,intercept=0, linetype="dashed")+
  theme(legend.position.inside = c(0.65,0.25),
        legend.background=element_rect(fill = alpha("white", 0))) +
  geom_text(AUC |> filter(toxinpresent=='Before'), mapping=aes(0.81,0.27, label=paste0('AUC = ',round(auc,2))), color='#88CCEE',size=3.4)+
  geom_text(AUC |> filter(toxinpresent=='During'), mapping=aes(0.81,0.19, label=paste0('AUC = ',round(auc,2))), color='#999933',size=3.4) +
  geom_text(AUC |> filter(toxinpresent=='After'), mapping=aes(0.81,0.12, label=paste0('AUC = ',round(auc,2))), color='#44AA99',size=3.4)


## plot confusion matrix ####
cm <- confusionMatrix(testing.dat$prediction, testing.dat$toxinpresent, dnn=c('Predicted', 'Observed'))

accuracy <- cm$overall

cm <- as.data.frame(cm$table) |>
  group_by(Observed) |>
  mutate(prop = Freq/sum(Freq),
         prop = round(prop,2))



cmplot<-ggplot(cm, aes(x = Observed, y = Predicted, fill=Freq)) +
  geom_tile(color="black") +
  scale_x_discrete(expand = c(0, 0))+ #remove white space
  scale_y_discrete(expand = c(0, 0))+ #remove white space
  scale_fill_gradient(low="white", high="#999933",
                      name="Frequency") +
  geom_text(aes(label = paste0("n=",Freq)), vjust = .5,  alpha = 1, size=3) +
  geom_text(aes(label = paste0("prop.=",prop)), vjust = 2.0,  alpha = 1, size=2) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y=element_text(angle=45),
        legend.position = 'none') +
  labs(title = paste0("Accuracy: ",round(accuracy,2))) 


cmplot+rocplot

ggsave('Figures/RandomForest/cm_roc_bda.png', height=4.5, width=6.5, dpi=1200)

# 6. look at top 12 variables ####


# top 4 have MDA>2 - present separately
top12 <- plt.dat |>
  slice(1:12)

bda_cyano <- sd_data  |>
  # we only want to keep data for RF that corresponds to the location and year a bloom occurred
  right_join(toxin_distance_mins |> dplyr::select(WaterbodyName, Year)) |>
  left_join(toxin_distance_mins) |>
  distinct() |>
  mutate(toxinpresent=ifelse(toxinpresent==1,'During',NA)) |>
  # i wish i could figure out how to code this instead of manually, but idk how
  mutate(toxinpresent=ifelse(WaterbodyName=='Cottonwood Creek Bay' &
                               is.na(toxinpresent),'Before',toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Fremont Bay' &
                               is.na(toxinpresent) &
                               month %in% c('May','Jun'),'Before', toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Fremont Bay' &
                               is.na(toxinpresent), 'After',toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Lacustrine Pelagic: Dam' &
                               is.na(toxinpresent) &
                               month %in% c('May','Jun', 'Jul'), 'Before', toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Lacustrine Pelagic: Dam' &
                               is.na(toxinpresent), 'After', toxinpresent),
         toxinpresent=ifelse(WaterbodyName=='Tough Creek Campground' &
                               is.na(toxinpresent), 'Before', toxinpresent)) |>
  # get rid of variables we don't need
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, -Cyanobacteria, -CHLA,-Latitude,-Longitude, - Lat, - Long, -ID, -distance_km) |>
  distinct() |>
  pivot_longer(cols=c(2:17), names_to = 'ShortName_Revised', values_to='ChemValue') |>
  filter(ShortName_Revised %in% c(top12$var)) |>
  mutate(ShortName_Revised = factor(ShortName_Revised, levels=c(top12$var))) |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  mutate(toxinpresent=factor(toxinpresent, levels=c('Before','During','After')))



library(ggpubr)
scales::viridis_pal(option='turbo')(7)
#"#30123BFF" "#4686FBFF" "#1AE4B6FF" "#A2FC3CFF" "#FABA39FF" "#E4460AFF" "#7A0403FF"

ggplot(bda_cyano, aes(toxinpresent, ChemValue)) +
  geom_boxplot() +
  geom_jitter(aes(color=WaterbodyName), alpha=0.5) +
  scale_color_viridis_d('',option='turbo') +
  stat_compare_means(fontface='bold',label = 'p.signif',comparisons = list(c('Before','During'), c('During','After'), c('Before','After'))) +   #  Kruskal-Wallis test 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # expands y-axis so we can see all results of kruskal-wallis comparisons
  facet_wrap(~ShortName_Revised, scales='free',ncol=4) +
  theme_bw() +
  labs(x='',y='')

## top 4 vars ####
bda_cyano |>
  filter(ShortName_Revised %in% c((top12 |> slice(1:4))$var)) |>
  data.frame() |>
  mutate(ShortName_Revised = factor(ShortName_Revised, levels=c(top12$var), labels=c(expression('SpC'~(µS~cm^-1)), expression('TN'~(mg~L^-1)), expression('Secchi depth'~(m)), expression('DO'~(mg~L^-1)), 'x','x','x','x','x','x','x','x'))) |>
  # this madness is just making pretty labels
  ggplot(aes(toxinpresent, ChemValue)) +
  geom_boxplot() +
  geom_jitter(aes(color=WaterbodyName), alpha=0.5) +
  scale_color_viridis_d('',option='turbo') +
  stat_compare_means(fontface='bold',label = 'p.signif',comparisons = list(c('Before','During'), c('During','After'), c('Before','After'))) +   #  Kruskal-Wallis test 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # expands y-axis so we can see all results of kruskal-wallis comparisons
  facet_wrap(~ShortName_Revised, scales='free',ncol=2,labeller=label_parsed) +
  theme_bw() +
  labs(x='',y='')
ggsave('Figures/RandomForest/top4boxplots.png', height=4.5, width=6.5, dpi=1200)

## next 8 vars ####
bda_cyano |>
  filter(ShortName_Revised %in% c((top12 |> slice(5:12))$var)) |>
  data.frame() |>
  mutate(ShortName_Revised = factor(ShortName_Revised, levels=c(top12$var), labels=c('x','x','x','x', expression('Phosphate'~(mg~L^-1)), expression('Inorganic N'~(mg~L^-1)), expression('Temp'~(degree*C)), 'IN:PO4~ratio', 'Phyto~Diversity~(H)', 'pH', expression('Schmidt Stability Index'~(J~m^-2)), 'TN:TP~ratio'))) |>
  # this madness is just making pretty labels
  ggplot(aes(toxinpresent, ChemValue)) +
  geom_boxplot() +
  geom_jitter(aes(color=WaterbodyName), alpha=0.5) +
  scale_color_viridis_d('',option='turbo') +
  stat_compare_means(fontface='bold',label = 'p.signif',comparisons = list(c('Before','During'), c('During','After'), c('Before','After'))) +   #  Kruskal-Wallis test 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # expands y-axis so we can see all results of kruskal-wallis comparisons
  facet_wrap(~ShortName_Revised, scales='free',ncol=4,labeller=label_parsed) +
  theme_bw() +
  labs(x='',y='')+
  theme(legend.position = 'bottom')
ggsave('Figures/RandomForest/last8boxplots.png', height=4.5, width=8.5, dpi=1200)


# 7. ID where/when are toxins occurring ####
cyano_checks <- cyano_prep |>
  st_as_sf(coords=c('Long','Lat'),crs=4326)

ggplot() +
  geom_sf(st_geometry(lake_shapefile), mapping=aes()) +
  geom_sf(cyano_checks |> filter(Year==2020), mapping=aes()) +
  geom_sf_text_repel(cyano_checks |> filter(Year==2020), mapping=aes(label=ID))
# 1 toxin producing bloom verified near LPD (2020-july)

ggplot() +
  geom_sf(st_geometry(lake_shapefile), mapping=aes()) +
  geom_sf(cyano_checks |> filter(Year==2021), mapping=aes()) +
  geom_sf_text_repel(cyano_checks |> filter(Year==2021), mapping=aes(label=ID))
# 7 toxin producing blooms verified. 3 months Aug-Oct near LPD; 4 months Jul-Oct near TCC

ggplot() +
  geom_sf(st_geometry(lake_shapefile), mapping=aes()) +
  geom_sf(cyano_checks |> filter(Year==2022), mapping=aes()) +
  geom_sf_text_repel(cyano_checks |> filter(Year==2022), mapping=aes(label=ID))
# 12 toxin producing blooms verified. by LPD, Aug-Oct (2 in sep); by TCC Aug-Oct (2 in Sep), by CCB Sep-Oct (2 ea.); by Fremont Bay 1 during July

ggplot() +
  geom_sf(st_geometry(lake_shapefile), mapping=aes()) +
  geom_sf(cyano_checks |> filter(Year==2023), mapping=aes()) +
  geom_sf_text_repel(cyano_checks |> filter(Year==2023), mapping=aes(label=ID))
# 8 toxin producing blooms verified. 3 months Aug-Oct near LPD (2 in sep, 2 in Oct, 1 aug); 3 months near CCB (1 ea.)  
