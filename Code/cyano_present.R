#------------------------------------------------#
# Script to look at conditions when toxin present 
#------------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. Quick look at parameters ####

nutrient_summarise <- BoysenNutrient |>
  group_by(CollDate, Year, ShortName_Revised) |>
  summarise(mean = mean(ChemValue),
            max = max(ChemValue),
            min = min(ChemValue)) |>
  ungroup()

chem_summarise <- BoysenChem |>
  group_by(CollDate, Year, ShortName_Revised) |>
  summarise(mean = mean(ChemValue),
            max = max(ChemValue),
            min = min(ChemValue)) |>
  ungroup()

data_summarise <- bind_rows(nutrient_summarise, chem_summarise) 


ggplot(data_summarise) +
  geom_point(aes(CollDate, mean)) +
  geom_errorbar(aes(CollDate, mean, ymin=min, ymax=max, width=0.2)) +
  geom_point(cyanotoxin |> filter(toxinpresent=='Y'), mapping=aes(CollDate, 0), shape=3, color='red4') +
  theme_bw() +
  labs(x='',y='Reservoir average and range') +
  facet_wrap(~ShortName_Revised, scales='free',nrow=5) 


# 2. Logistic regression for cyanos present/not present ####
cyano_prep <- cyanotoxin |>
  mutate(month = month(CollDate,label=TRUE,abbr=TRUE),
         Year=year(CollDate)) |>
  select(-CollDate, -Advisory) |>
  mutate(toxinpresent=as.character(toxinpresent)) |>
  mutate(toxinpresent=ifelse(toxinpresent=='N',0,1)) |>
  filter(toxinpresent != 0) |>
  distinct()



logistic_data <- data_summarise |>
   mutate(month = month(CollDate,label=TRUE,abbr=TRUE),
          Year=year(CollDate)) |>
  select(-CollDate, -max, -min) |>
  pivot_wider(names_from = ShortName_Revised, values_from = mean) |>
  left_join(cyano_prep) |>
  mutate(toxinpresent=ifelse(is.na(toxinpresent), 0, toxinpresent))




# m1 <- glm(toxinpresent~CHLA, logistic_data, family=binomial)
# summary(m1)

m1 <- glm(toxinpresent~DO, logistic_data, family=binomial)
summary(m1)

ggplot(logistic_data, mapping = aes(DO, toxinpresent)) +
  geom_point(shape = 1) +
  stat_smooth(method = 'glm', se = FALSE, fullrange = TRUE, method.args = list(family = binomial)) +
  theme_bw() +
  labs(x = 'DO'~(mg~L^-1),
       y = 'Cyanotoxins present')

# m2 <- glm(toxinpresent~H, logistic_data, family=binomial)
# summary(m2)

# m2 <- glm(toxinpresent~IN.PO4, logistic_data, family=binomial)
# summary(m2)

# m2 <- glm(toxinpresent~NH4, logistic_data, family=binomial)
# summary(m2)

# m2 <- glm(toxinpresent~NO3, logistic_data, family=binomial)
# summary(m2)

m2 <- glm(toxinpresent~pH, logistic_data, family=binomial)
summary(m2)

ggplot(logistic_data, mapping = aes(pH, toxinpresent)) +
  geom_point(shape = 1) +
  stat_smooth(method = 'glm', se = FALSE, fullrange = TRUE, method.args = list(family = binomial)) +
  theme_bw() +
  labs(x = 'pH',
       y = 'Cyanotoxins present')


# m3 <- glm(toxinpresent~PO4, logistic_data, family=binomial)
# summary(m3)

m3 <- glm(toxinpresent~Secchi, logistic_data, family=binomial)
summary(m3)

ggplot(logistic_data, mapping = aes(Secchi, toxinpresent)) +
  geom_point(shape = 1) +
  stat_smooth(method = 'glm', se = FALSE, fullrange = TRUE, method.args = list(family = binomial)) +
  theme_bw() +
  labs(x = 'Secchi depth (m)',
       y = 'Cyanotoxins present')


# m4 <- glm(toxinpresent~SpC, logistic_data, family=binomial)
# summary(m4)

# m4 <- glm(toxinpresent~Stability, logistic_data, family=binomial)
# summary(m4)

# m4 <- glm(toxinpresent~Temp, logistic_data, family=binomial)
# summary(m4)

# m4 <- glm(toxinpresent~TN, logistic_data, family=binomial)
# summary(m4)

# m4 <- glm(toxinpresent~TN.TP, logistic_data, family=binomial)
# summary(m4)

# m4 <- glm(toxinpresent~TP, logistic_data, family=binomial)
# summary(m4)

# m.4 <- glm(toxinpresent~DO*Secchi, logistic_data, family=binomial)
# summary(m.4)

# m.4 <- glm(toxinpresent~DO*pH, logistic_data, family=binomial)
# summary(m.4)


# I don't want to lose everything by dropping NAs, but the stepwise model finder can't handle NAs :( 

# logistic_data_nona <- drop_na(logistic_data)
# # fit intercept model
# intercept.model<-glm(toxinpresent ~ 1, data=logistic_data_nona, family=binomial)
# 
# # perform forward selection using BIC (with k=log(n))
# step.modelBIC<-MASS::stepAIC(intercept.model, direction="forward",
#                              scope = ~ CHLA + DO,
#                              trace=F, k=log(27))
# 
# # perform forward selection using BIC (with k=log(n))
# step.modelBIC<-MASS::stepAIC(intercept.model, direction="forward",
#                              scope = ~ CHLA + DO + H + IN.PO4 + NH4 + NO3 + pH + PO4 + Secchi + SpC + Stability + Temp + TN + TN.TP + TP,
#                              trace=F, k=log(27))
# 
# summary(step.modelBIC)
# 
# 
# # display the summary table for the model chosen by forward BIC selection
# broom::tidy(step.modelBIC) %>%
#   mutate(p.value = scales::pvalue(p.value, accuracy = 0.0001)) %>% knitr::kable(
#     caption = "Estimates for the MLR model selected 
#                 through forward selection using BIC as criterion.
#                 The model is showing the total number of daily bike rentals 
#                 as a function of the listed predictors.",
#     col.names = c("Predictor", "B", "SE", "t", "p"),
#     digits = c(3, 3, 3, 3, 4),
#     align = c("l", "r", "r", "r", "r"),
#     format.args = list(big.mark = ",")
#   )


# 3. Before, during, after Cyanos ####
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
  mutate(toxinpresent=factor(toxinpresent, levels=c('Before','During','After')))

library(ggpubr)
ggplot(bda_cyano) +
  geom_boxplot(aes(toxinpresent, ChemValue)) +
  geom_jitter(aes(toxinpresent, ChemValue, color=WaterbodyName), alpha=0.5) +
  scale_color_viridis_d(option='turbo') +
  stat_compare_means(fontface='bold',label = 'p.signif',comparisons = list(c('Before','During'), c('During','After'), c('Before','After'))) +   #  Kruskal-Wallis test 
  facet_wrap(~ShortName_Revised, scales='free')
  
  






