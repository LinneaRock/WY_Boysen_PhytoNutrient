#--------------------------------------------------------------------------#
# Structural equation modeling to predict high percentages of cyanobacteria
#--------------------------------------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. load libraries needed specifically for this script ####
library(tidySEM)
library(lavaan)
library(statpsych) # for skew.test

# 2. Prep data ####
# create metadata 
sd_data <- BoysenNutrient |>
  bind_rows(BoysenChem |> select(-maxdepth)) |>
  left_join(BoysenChem |> select(WaterbodyName, CollDate, maxdepth) |> distinct()) |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate, maxdepth) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) 

# look at correlations 
psych::pairs.panels(sd_data |> select(-c(1:8)))

# 3. Normality assumptions test ####
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

## 3a. transform to normality ####
# for biomass % and all other non-normal variables (p<0.05 in skew dataset), add1 and log10

trn_data <- sd_data |>
  mutate_at(vars(Diatom, `Green algae`, Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate, CHLA, TN.TP, Stability, IN.PO4, NO3, DO, TN, maxdepth), ~log10(.+1))
  # scale the data
 


# 4. SEM ####

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
