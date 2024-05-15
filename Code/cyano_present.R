#------------------------------------------------#
# Script to look at conditions when toxin present 
#------------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. Quick look at parameters ####

nutrient_summarise <- BoysenNutrient |>
  select(-ChemUnits, -SampleDepth) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) |>
  mutate(IN=NO3+NH4) |>
  pivot_longer(c(TN.TP, IN.PO4, TN, TP, PO4, NO3, NH4, IN, CHLA), names_to = 'ShortName_Revised', values_to = 'ChemValue') |>
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

nutrient_summarise |>
  filter(ShortName_Revised %in% c('IN', 'TN')) |>
  left_join(cyanotoxin) |>
ggplot() +
  geom_point(aes(CollDate, mean)) +
  geom_errorbar(aes(CollDate, mean, ymin=min, ymax=max, width=0.2)) +
  facet_wrap(~ShortName_Revised, scales='free', labeller=as_labeller(c(IN='Inorganic N', TN='Total N'))) +
  theme_bw() +
  annotate("rect", xmin = as.Date('2020-01-01'), xmax = as.Date('2020-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
  annotate("rect", xmin = as.Date('2020-10-15'), xmax = as.Date('2021-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
  annotate("rect", xmin = as.Date('2021-10-15'), xmax = as.Date('2022-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
  annotate("rect", xmin = as.Date('2022-10-15'), xmax = as.Date('2023-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
 # scale_color_manual('', values=c('#882255','#332288')) +
  theme(legend.position = 'none') +
  labs(x='',y='Concentration'~(mg~L^-1))
  ggsave('Figures/in-res_N.png',height=4.5,width=6.5, units='in',dpi=1200)
  
  
  
  nutrient_summarise |>
    filter(ShortName_Revised %in% c('PO4', 'TP')) |>
    ggplot() +
    geom_point(aes(CollDate, mean)) +
    geom_errorbar(aes(CollDate, mean, ymin=min, ymax=max, width=0.2)) +
    facet_wrap(~ShortName_Revised, scales='free', labeller=as_labeller(c(PO4='Inorganic P', TP='Total P'))) +
    theme_bw() +
    annotate("rect", xmin = as.Date('2020-01-01'), xmax = as.Date('2020-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
    annotate("rect", xmin = as.Date('2020-10-15'), xmax = as.Date('2021-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
    annotate("rect", xmin = as.Date('2021-10-15'), xmax = as.Date('2022-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
    annotate("rect", xmin = as.Date('2022-10-15'), xmax = as.Date('2023-04-15'), ymin = -Inf, ymax = Inf, alpha = 0.5, color = "grey") +
   # scale_color_manual('', values=c('#117733','#DDCC77')) +
    theme(legend.position = 'none') +
    labs(x='',y='Concentration'~(mg~L^-1))
  ggsave('Figures/in-res_P.png',height=4.5,width=6.5, units='in',dpi=1200)




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

 # m.4 <- glm(toxinpresent~maxdepth, logistic_data, family=binomial)
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
ggplot(bda_cyano, aes(toxinpresent, ChemValue)) +
  geom_boxplot() +
  geom_jitter(aes(color=WaterbodyName), alpha=0.5) +
  scale_color_viridis_d(option='turbo') +
  stat_compare_means(fontface='bold',label = 'p.signif',comparisons = list(c('Before','During'), c('During','After'), c('Before','After'))) +   #  Kruskal-Wallis test 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # expands y-axis so we can see all results of kruskal-wallis comparisons
  facet_wrap(~ShortName_Revised, scales='free')
  






