#--------------------------------------------------------------------------#
# Various models to predict high percentages of cyanobacteria
#--------------------------------------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. Structural Equation Modeling ####

## 1A. load libraries needed specifically for this script ####
library(tidySEM)
library(lavaan)
library(statpsych) # for skew.test

## 1B. Prep data ####
sd_data <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) 

# look at correlations 
psych::pairs.panels(sd_data |> select(-c(1:8)))

## 1C. Normality assumptions test ####
variables <- colnames(sd_data |> select(-c(1:8)))

skew <- data.frame(Skewness=NA, p=NA, variable=NA)

for(var in variables) {
  dat <- sd_data |> select(var) |>
    drop_na() 
  
  dat <- as.vector(dat |> select(1) |> distinct())[[1]]
  
  tmp <- test.skew(dat) |>
    as.data.frame() |>
    mutate(variable=var)
  
  skew <- bind_rows(skew, tmp)
}

## 1D. transform to normality ####
# for biomass % and all other non-normal variables (p<0.05 in skew dataset), add1 and log10

trn_data <- sd_data |>
  mutate_at(vars(Diatom, `Green algae`, Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate, CHLA, TN.TP, Stability, IN.PO4, NO3, DO, TN, maxdepth), ~log10(.+1))
  # scale the data
 


# 1E. Run the SEM ####

# starting model based on Deutsch et al., 2020
m1 <- 'Cyanobacteria ~ Temp + Stability + maxdepth + TP + TN + TN.TP + IN.PO4 + CHLA
Stability ~ Temp + DO + maxdepth
CHLA ~ Temp + Secchi + Stability + pH
DO ~ CHLA'

fit1 <- sem(m1, trn_data)
summary(fit1, standardized = TRUE)
graph_sem(fit1)
semPlot::semPaths(fit1, what = "std", whatLabels = "std", residuals=FALSE)
modificationindices(fit1)





# testing stratification model
m2 <- 'Stratification=~Secchi + Stability + Temp
       Secchi~~maxdepth
       Stability~Temp + maxdepth
       Cyanobacteria ~ Stratification'

fit2 <- sem(m2, trn_data)
summary(fit2, standardized = TRUE)
summary(fit2, modindices=TRUE)
graph_sem(fit2)
semPlot::semPaths(fit2, what = "std", whatLabels = "std", residuals=FALSE)
modificationindices(fit2)





# 2. CART modeling ####

## 2A. load packages ####

library(tidyverse)
library(rpart)
library(rsample)
library(caret)
library(rpart.plot)
library(yardstick)



## 2B. Predict Before-During-After toxins present ####
cart_data <- trn_data  |>
  left_join(cyano_prep) |> # this df is from cyano_present.R
  mutate(toxinpresent=ifelse(toxinpresent==1,'During',NA)) |>
    # i wish i could figure out how to code this instead of manually, but idk how
  # 2020
  mutate(toxinpresent=if_else(Year==2020 & month%in%c('May','Jun'),'Before',
                              ifelse(Year==2020 & is.na(toxinpresent), 'After', toxinpresent))) |>
  # 2021
  mutate(toxinpresent=if_else(Year==2021 & month%in%c('May','Jun'),'Before',toxinpresent)) |>
  # 2022
  mutate(toxinpresent=if_else(Year==2022 & month%in%c('May','Jun'),'Before',toxinpresent)) |>
  # 2023
  mutate(toxinpresent=if_else(Year==2023 & month%in%c('May','Jun','Jul'),'Before',toxinpresent)) |>
  mutate(toxinpresent=factor(toxinpresent, levels=c('Before','During','After')))|>
  # time is too important, so get rid of it to assess other variables
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, -Cyanobacteria)



# Training data 
# split training and testing data by 60/40% 
test_inds <- initial_split(cart_data, 0.6)
# Split data into test/train 
train <- training(test_inds) 
test <- testing(test_inds) 


# Run CART model
CART_mod <- rpart(toxinpresent~., data=train, method='class', cp=0.01)

# test cart model 
test$predict <- predict(CART_mod, test, 'class')
cm <- conf_mat(test, toxinpresent, predict)
accuracy(test, toxinpresent, predict) # like r2
cm

## get the best cp 
printcp(CART_mod)
bestcp <- CART_mod$cptable[which.min(CART_mod$cptable[,'xerror']), 'CP']

CART_mod_pruned <- prune(CART_mod, cp=bestcp)

## print tree 

plot(CART_mod_pruned)
text(CART_mod_pruned, cex=0.8, xpd=TRUE)
#text(CART_mod_pruned, cex=0.8, use.n=TRUE, xpd=TRUE)



## 2C. Predict toxin presence/absence ####
cart_data <- trn_data  |>
  left_join(cyano_prep) |>
  mutate(toxinpresent = ifelse(is.na(toxinpresent), 0, toxinpresent)) |>
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, -Cyanobacteria) |>
  mutate(toxinpresent=as.factor(toxinpresent))

# Training data 
# split training and testing data by 60/40% 
test_inds <- initial_split(cart_data, 0.6)
# Split data into test/train 
train <- training(test_inds) 
test <- testing(test_inds) 



# Run CART model 
CART_mod <- rpart(toxinpresent~., data=train, method='class', cp=0.01)

# test cart model
test$predict <- predict(CART_mod, test, 'class')
cm <- conf_mat(test, toxinpresent, predict)
accuracy(test, toxinpresent, predict) # like r2
cm

## get the best cp 
printcp(CART_mod)
bestcp <- CART_mod$cptable[which.min(CART_mod$cptable[,'xerror']), 'CP']

CART_mod_pruned <- prune(CART_mod, cp=bestcp)

## print tree 

plot(CART_mod_pruned)
text(CART_mod_pruned, cex=0.8, xpd=TRUE)
#text(CART_mod_pruned, cex=0.8, use.n=TRUE, xpd=TRUE)




## 2D. Predict high vs low cyanobacteria (greater than 50% biomass) ####
cart_data <- trn_data |>
  mutate(cyano= ifelse(Cyanobacteria > log10(51), 'high cyano','low cyano'))|>
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, -Cyanobacteria) |>
  mutate(cyano=as.factor(cyano))

# Training data 
# split training and testing data by 60/40% 
test_inds <- initial_split(cart_data, 0.6)
# Split data into test/train 
train <- training(test_inds) 
test <- testing(test_inds) 



# Run CART model 
CART_mod <- rpart(cyano~., data=train, method='class', cp=0.01)

## test cart model 
test$predict <- predict(CART_mod, test, 'class')
cm <- conf_mat(test, cyano, predict)
accuracy(test, cyano, predict) # like r2
cm

## get the best cp 
printcp(CART_mod)
bestcp <- CART_mod$cptable[which.min(CART_mod$cptable[,'xerror']), 'CP']

CART_mod_pruned <- prune(CART_mod, cp=bestcp)

## print tree 

plot(CART_mod_pruned)
text(CART_mod_pruned, cex=0.8, xpd=TRUE)
#text(CART_mod_pruned, cex=0.8, use.n=TRUE, xpd=TRUE)



# 3. Random Forest modeling ####

## 3A. load packages ####
library(randomForest)
library(rsample)
library(caret)
library(reprtree)


# RF can't handle NAs, so let's remove the CHLA, it hasn't been a significant predictor of anything yet anyway

# data do not need to be normal for RF - and scaling/transforming data can be harmful


## 3B. Regression RF to predict cyanobacteria % biomass ####

rf.data <- sd_data |>
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, - CHLA) 

# use 80% data for training, 20% for testing
set.seed(62693)
split <- initial_split(rf.data, prop=0.80)
training.dat <- training(split) |> mutate_if(is.numeric, round, digits=2)
testing.dat <- testing(split) |> mutate_if(is.numeric, round, digits=2) 


# model with default parameters
rf_default <- train(Cyanobacteria ~.,
                    training.dat,
                    metric='RMSE',
                    method='rf',
                    tuneGrid=expand.grid(.mtry=ncol(training.dat)/3),
                    ntree=500,
                    trControl=trainControl(method='cv', number=10))

rf_default #rmse 25.12346



# find best mtry
set.seed(62693)
rf_mtry <- train(Cyanobacteria~.,
                 data = training.dat,
                 method = "rf",
                 metric = "RMSE",
                 tuneGrid = expand.grid(.mtry = c(1: 10)),
                 trControl = trainControl(method = "cv",
                                          number = 10,
                                          search = "grid"))

print(rf_mtry) 
plot(rf_mtry) 
# mtry=8, with RMSE of  25.06891 


# find best ntrees
store_maxtrees <- list()
for (ntree in c(100, 150, 250, 300, 350, 400, 450, 500, 800, 1000, 2000)) {
  set.seed(62693)
  rf_maxtrees <- train(Cyanobacteria~.,
                       data = training.dat,
                       method = "rf",
                       metric = "RMSE",
                       tuneGrid = expand.grid(.mtry = c(1: 10)),
                       trControl = trainControl(method = "cv",
                                                number = 10,
                                                search = "grid"),
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)

summary(results_tree) # 100



set.seed(06261993)
# fit the model with the best hyperparameters
BIOMASS_fit_rf <- randomForest(Cyanobacteria~.,
                       training.dat,
                       method = "rf",
                       metric = "RMSE",
                       tuneGrid = expand.grid(.mtry = c(1: 10)),
                       trControl = trainControl(method = "cv",
                                                number = 10),
                       importance = TRUE,
                       mtry = 8,
                       ntree = 100)
# get predicted values
testing.dat$prediction <- predict(BIOMASS_fit_rf, testing.dat)



BIOMASS_fit_rf

# Call:
#   randomForest(formula = Cyanobacteria ~ ., data = training.dat,      method = "rf", metric = "RMSE", tuneGrid = expand.grid(.mtry = c(1:10)),      trControl = trainControl(method = "cv", number = 10), importance = TRUE,      mtry = 8, ntree = 100) 
# Type of random forest: regression
# Number of trees: 100
# No. of variables tried at each split: 8
# 
# Mean of squared residuals: 621.9576
# % Var explained: 59.4

varImpPlot(BIOMASS_fit_rf)
plot(BIOMASS_fit_rf)

plt.dat<-varImpPlot(BIOMASS_fit_rf, type = 1, scale = TRUE,
           n.var = ncol(rf.data) - 1, cex = 0.8,
           main = "Variable importance") |>
  as.data.frame() |>
  arrange(desc(`%IncMSE`))

plt.dat<-cbind(rownames(plt.dat), plt.dat) |>
  rename(var=1)

ggplot(plt.dat |> filter(var != 'NH4')) +
  geom_bar(aes(x=`%IncMSE`,y=reorder(var, `%IncMSE`)), stat='identity') +
  theme_minimal() +
  labs(y='')

ggsave('Figures/RandomForest/Cyano_%biomass/Var_importance_biomass.png', height=4.5, width=6.5, dpi=1200)

BIOMASS_fit_rf$importance

library(reprtree)
plot.getTree(BIOMASS_fit_rf)



BIOMASS_fit_rf$mse[length(BIOMASS_fit_rf$mse)]
# take square root to calculate RMSE for the model
sqrt(BIOMASS_fit_rf$mse[length(BIOMASS_fit_rf$mse)]) # 24.93908

# now illustrate how to calculate RMSE on test data vs. training data
predValues <- predict(BIOMASS_fit_rf,testing.dat)
# we can calculate it  directly 
sqrt(mean((testing.dat$Cyanobacteria - predValues)^2)) # 18.74332



test.result<- lm(Cyanobacteria~prediction, testing.dat)
summary(test.result) # Adjusted R-squared:  0.6834


ggplot() +
  geom_point(testing.dat, mapping=aes(prediction, Cyanobacteria)) +
  theme_classic() +
  labs(x='Predicted % biomass cyanobacteria', y='Observed % biomass cyanobacteria') +
  geom_abline(slope=1, intercept=0, color='red4') +
  geom_text(mapping=aes(25,95, label=paste0('y = ', round(coef(test.result)[[2]], 2), 'x ', round(coef(test.result)[[1]], 2)))) +
  geom_text(mapping=aes(25,90, label=paste0('R² = ', round(summary(test.result)$adj.r,2)))) +
  geom_text(mapping=aes(25,85, label=paste0('RMSE = ', round(sqrt(mean((testing.dat$Cyanobacteria - predValues)^2))))))

ggsave('Figures/RandomForest/Cyano_%biomass/prediction_line.png', height=4.5, width=6.5, dpi=1200)






## 3C. Categorical RF to predict presence absence ####
rf.data <- sd_data  |>
  left_join(cyano_prep) |>
  mutate(toxinpresent = ifelse(is.na(toxinpresent), 0, toxinpresent)) |>
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, -Cyanobacteria,-CHLA) |>
  mutate(toxinpresent=as.factor(toxinpresent))

# use 80% data for training, 20% for testing
set.seed(62693)
split <- initial_split(rf.data, prop=0.80)
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

rf_default #accuracy=0.8857



# find best mtry
set.seed(62693)
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
# mtry=6, with RMSE of 0.8847


# find best ntrees
store_maxtrees <- list()
for (ntree in c(100, 150, 250, 300, 350, 400, 450, 500, 800, 1000, 2000)) {
  set.seed(62693)
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

summary(results_tree) # they're like all the same? so I'll just go with 100



set.seed(06261993)
# fit the model with the best hyperparameters
presence_fit_rf <- randomForest(toxinpresent~.,
                       training.dat,
                       method = "rf",
                       #metric = "RMSE",
                       tuneGrid = expand.grid(.mtry = c(1: 10)),
                       trControl = trainControl(method = "cv",
                                                number = 10),
                       importance = TRUE,
                       mtry = 6,
                       ntree = 100)
# get predicted values
testing.dat$prediction <- predict(presence_fit_rf, testing.dat)



presence_fit_rf

# Call:
#   randomForest(formula = toxinpresent ~ ., data = training.dat,      method = "rf", tuneGrid = expand.grid(.mtry = c(1:10)), trControl = trainControl(method = "cv",          number = 10), importance = TRUE, mtry = 6, ntree = 100) 
# Type of random forest: classification
# Number of trees: 100
# No. of variables tried at each split: 6
# 
# OOB estimate of  error rate: 11.45%
# Confusion matrix:
#    0   1 class.error
# 0  58  7   0.1076923
# 1  8  58   0.1212121

varImpPlot(presence_fit_rf)
plot(presence_fit_rf)

varImpPlot(presence_fit_rf, type = 1, scale = TRUE,
           n.var = ncol(rf.data) - 1, cex = 0.8,
           main = "Variable importance")

presence_fit_rf$importance

plt.dat<-varImpPlot(presence_fit_rf, type = 1, scale = TRUE,
                    n.var = ncol(rf.data) - 1, cex = 0.8,
                    main = "Variable importance") |>
  as.data.frame() |>
  arrange(desc(MeanDecreaseAccuracy))

plt.dat<-cbind(rownames(plt.dat), plt.dat) |>
  rename(var=1)

ggplot(plt.dat |> filter(var != 'NH4')) +
  geom_bar(aes(x=MeanDecreaseAccuracy,y=reorder(var, MeanDecreaseAccuracy)), stat='identity') +
  theme_minimal() +
  labs(y='',x='Variable importance')

ggsave('Figures/RandomForest/PresenceToxin/Var_importance_presence.png', height=4.5, width=6.5, dpi=1200)



library(reprtree)
plot.getTree(presence_fit_rf)


# # plot ROC curve
# # predict test set, get probs instead of response
# predictions <- as.data.frame(predict(presence_fit_rf, testing.dat, type = "prob")) 
# 
# 
# # predict class and then attach test class
# predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
# predictions$observed <- testing.dat$toxinpresent
# head(predictions)
# 
# library(pROC)
# roc.absent <- roc(ifelse(predictions$observed==0,'absent',"non-absent"), as.numeric(predictions$`0`))
# 
# roc.present <- roc(ifelse(predictions$observed==1, "present", "non-present"), as.numeric(predictions$`0`))
# 
# 
# 
# roc.dat <- data.frame(toxinpresent=0,
#                       Sensitivity=roc.absent[["sensitivities"]],
#                       Specificity=roc.absent[["specificities"]]) |>
#   rbind(data.frame(toxinpresent=1,
#                    Sensitivity=roc.present[["sensitivities"]],
#                    Specificity=roc.present[["specificities"]])) 
# 
# AUC <- data.frame(toxinpresent=c('Before','During','After'),
#                   auc=c(roc.Before[["auc"]],roc.During[["auc"]],roc.After[["auc"]]))
# 
# 
# 
# ggplot() +
#   theme_classic() +
#   geom_path(roc.dat, mapping=aes(x=1-Specificity,y=Sensitivity,color = toxinpresent)) +
#   scale_color_manual('Cyanotoxin presence', values=c('#88CCEE','#999933')) +
#   geom_abline(slope=1,intercept=0, linetype="dashed")+
#   theme(legend.position = c(0.65,0.25),
#         legend.background=element_rect(fill = alpha("white", 0))) +
#   geom_text(AUC |> filter(toxinpresent=='Before'), mapping=aes(0.81,0.27, label=paste0('AUC = ',round(auc,2))), color='#88CCEE',size=3.4)+
#   geom_text(AUC |> filter(toxinpresent=='During'), mapping=aes(0.81,0.19, label=paste0('AUC = ',round(auc,2))), color='#999933',size=3.4) +
#   geom_text(AUC |> filter(toxinpresent=='After'), mapping=aes(0.81,0.12, label=paste0('AUC = ',round(auc,2))), color='#44AA99',size=3.4)
# 
# 
# 
# 
# # plot confusion matrix
# cm <- confusionMatrix(testing.dat$prediction, testing.dat$toxinpresent, dnn=c('Predicted', 'Observed'))
# 
# accuracy <- cm$overall
# 
# cm <- as.data.frame(cm$table) |>
#   group_by(Observed) |>
#   mutate(prop = Freq/sum(Freq),
#          prop = round(prop,2))
# 
# 
# 
# cmplot<-ggplot(cm, aes(x = Observed, y = Predicted, fill=Freq)) +
#   geom_tile(color="black") +
#   scale_x_discrete(expand = c(0, 0))+ #remove white space
#   scale_y_discrete(expand = c(0, 0))+ #remove white space
#   scale_fill_gradient(low="white", high="#999933",
#                       name="Frequency") +
#   geom_text(aes(label = paste0("n=",Freq)), vjust = .5,  alpha = 1, size=3) +
#   geom_text(aes(label = paste0("prop.=",prop)), vjust = 2.0,  alpha = 1, size=2) +
#   theme_minimal() +
#   theme(axis.title = element_blank(),
#         axis.text.y=element_text(angle=45),
#         legend.position = 'none') +
#   labs(title = paste0("Accuracy: ",round(accuracy,2))) 
# 
# library(patchwork)
# cmplot+rocplot
# 
# ggsave('Figures/RandomForest/beforeduringaftertoxin/cm_roc_bda.png', height=4.5, width=6.5, dpi=1200)












## 3D. Categorical RF to predict before-during-after toxin ####

rf.data <- sd_data  |>
  left_join(cyano_prep) |>
  mutate(toxinpresent=ifelse(toxinpresent==1,'During',NA)) |>
  # i wish i could figure out how to code this instead of manually, but idk how
  # 2020
  mutate(toxinpresent=if_else(Year==2020 & month%in%c('May','Jun'),'Before',
                              ifelse(Year==2020 & is.na(toxinpresent), 'After', toxinpresent))) |>
  # 2021
  mutate(toxinpresent=if_else(Year==2021 & month%in%c('May','Jun'),'Before',toxinpresent)) |>
  # 2022
  mutate(toxinpresent=if_else(Year==2022 & month%in%c('May','Jun'),'Before',toxinpresent)) |>
  # 2023
  mutate(toxinpresent=if_else(Year==2023 & month%in%c('May','Jun','Jul'),'Before',toxinpresent)) |>
  mutate(toxinpresent=factor(toxinpresent, levels=c('Before','During','After')))|>
  # time is too important, so get rid of it to assess other variables
  select(-Group, -CollDate,-Year,-month,-julianday,-Diatom,-`Green algae`,-Dinoflagellate, -`Golden algae`, -Flagellate, -Cyanobacteria, -CHLA)





# use 80% data for training, 20% for testing
set.seed(62693)
split <- initial_split(rf.data, prop=0.80)
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

rf_default #accuracy=0.8796



# find best mtry
set.seed(62693)
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
# mtry=4, with RMSE of 0.9010


# find best ntrees
store_maxtrees <- list()
for (ntree in c(100, 150, 250, 300, 350, 400, 450, 500, 800, 1000, 2000)) {
  set.seed(62693)
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

summary(results_tree) # 250



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
                       mtry = 4,
                       ntree = 250)
# get predicted values
testing.dat$prediction <- predict(bda_fit_rf, testing.dat)



bda_fit_rf

# Call:
#   randomForest(formula = toxinpresent ~ ., data = training.dat,      method = "rf", tuneGrid = expand.grid(.mtry = c(1:10)), trControl = trainControl(method = "cv",          number = 10), importance = TRUE, mtry = 4, ntree = 250) 
# Type of random forest: classification
# Number of trees: 250
# No. of variables tried at each split: 4
# 
# OOB estimate of  error rate: 12.21%
# Confusion matrix:
#         Before During After class.error
# Before     41      7     1  0.16326531
# During      6     60     0  0.09090909
# After       0      2    14  0.12500000

varImpPlot(bda_fit_rf)
plot(bda_fit_rf)

varImpPlot(bda_fit_rf, type = 1, scale = TRUE,
           n.var = ncol(rf.data) - 1, cex = 0.8,
           main = "Variable importance")

bda_fit_rf$importance


plot.getTree(bda_fit_rf)


plt.dat<-varImpPlot(bda_fit_rf, type = 1, scale = TRUE,
                    n.var = ncol(rf.data) - 1, cex = 0.8,
                    main = "Variable importance") |>
  as.data.frame() |>
  arrange(desc(MeanDecreaseAccuracy))

plt.dat<-cbind(rownames(plt.dat), plt.dat) |>
  rename(var=1)

ggplot(plt.dat) +
  geom_bar(aes(x=MeanDecreaseAccuracy,y=reorder(var, MeanDecreaseAccuracy)), stat='identity') +
  theme_minimal() +
  labs(y='',x='Mean decrease accuracy', title='Variable Importance')

ggsave('Figures/RandomForest/beforeduringaftertoxin/Var_importance_bda.png', height=4.5, width=6.5, dpi=1200)




# plot ROC curve
# predict test set, get probs instead of response
predictions <- as.data.frame(predict(bda_fit_rf, testing.dat, type = "prob"))

# predict class and then attach test class
predictions$predict <- names(predictions)[1:3][apply(predictions[,1:3], 1, which.max)]
predictions$observed <- testing.dat$toxinpresent
head(predictions)


library(pROC)
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
  theme(legend.position = c(0.65,0.25),
        legend.background=element_rect(fill = alpha("white", 0))) +
  geom_text(AUC |> filter(toxinpresent=='Before'), mapping=aes(0.81,0.27, label=paste0('AUC = ',round(auc,2))), color='#88CCEE',size=3.4)+
  geom_text(AUC |> filter(toxinpresent=='During'), mapping=aes(0.81,0.19, label=paste0('AUC = ',round(auc,2))), color='#999933',size=3.4) +
  geom_text(AUC |> filter(toxinpresent=='After'), mapping=aes(0.81,0.12, label=paste0('AUC = ',round(auc,2))), color='#44AA99',size=3.4)
  



# plot confusion matrix
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

library(patchwork)
cmplot+rocplot

ggsave('Figures/RandomForest/beforeduringaftertoxin/cm_roc_bda.png', height=4.5, width=6.5, dpi=1200)


### 3D-gpt ideas #### - pdp not working 
# librayr(pdp)
# # Partial Dependence Plot for specific conductivity
# pdp_spc<- partial(bda_fit_rf, pred.var = "SpC", plot = TRUE, rug = TRUE)
# # PDP for specific conductivity
# pdp_conductivity <- partial(bda_fit_rf, pred.var = "SpC", train = training.dat)
# p1 <- autoplot(pdp_conductivity) +
#   labs(title = "Partial Dependence of Specific Conductivity", x = "Specific Conductivity", y = "Predicted Response")
# 
# 
# 
# # Partial Dependence Plot for pH
# pdp_ph <- partial(rf_model, pred.var = "pH", plot = TRUE, rug = TRUE)
# # Partial Dependence Plot for water temperature
# pdp_temp <- partial(rf_model, pred.var = "water_temperature", plot = TRUE, rug = TRUE)




### 3D - variables ####
# look at top ten most important variables via Mean Decrease in Accuracy (MDA)
top10 <- plt.dat |>
  slice(1:10)

bda_cyano <- rbind(BoysenChem, BoysenNutrient) |>
  left_join(cyano_prep) |>
  mutate(toxinpresent=ifelse(toxinpresent==1,'During',NA)) |>
  # i wish i could figure out how to code this instead of manually, but idk how
  # 2020
  mutate(toxinpresent=if_else(Year==2020 & month%in%c('May','Jun'),'Before',
                              ifelse(Year==2020 & is.na(toxinpresent), 'After', toxinpresent))) |>
  # 2021
  mutate(toxinpresent=if_else(Year==2021 & month%in%c('May','Jun'),'Before',toxinpresent)) |>
  # 2022
  mutate(toxinpresent=if_else(Year==2022 & month%in%c('May','Jun'),'Before',toxinpresent)) |>
  # 2023
  mutate(toxinpresent=if_else(Year==2023 & month%in%c('May','Jun','Jul'),'Before',toxinpresent)) |>
  mutate(toxinpresent=factor(toxinpresent, levels=c('Before','During','After'))) |>
  # here, let's just keep the top 10 best predictor variables from MDA
  filter(ShortName_Revised %in% c(top10$var)) |>
  mutate(ShortName_Revised = factor(ShortName_Revised, levels=c(top10$var)))

#library(ggpubr)
ggplot(bda_cyano, aes(toxinpresent, ChemValue)) +
  geom_boxplot() +
  geom_jitter(aes(color=WaterbodyName), alpha=0.5) +
  scale_color_viridis_d('',option='turbo') +
  # stat_compare_means(fontface='bold',label = 'p.signif',comparisons = list(c('Before','During'), c('During','After'), c('Before','After'))) +   #  Kruskal-Wallis test 
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # expands y-axis so we can see all results of kruskal-wallis comparisons
  facet_wrap(~ShortName_Revised, scales='free',ncol=3) +
  theme_bw() +
  labs(x='',y='')
  

