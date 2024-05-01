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
  bind_rows(BoysenChem) |>
  left_join(BoysenPhyto_cat) |>
  mutate(Group = paste(WaterbodyName, CollDate, sep=' ')) |>
  select(Group, WaterbodyName, CollDate, Year, month, julianday, Latitude, Longitude, ShortName_Revised, ChemValue, Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate) |>
  pivot_wider(names_from=ShortName_Revised, values_from=ChemValue) 

# look at correlations 
psych::pairs.panels(sd_data |> select(-c(1:8)))

# 3. Normality assumptions test ####
variables <- colnames(sd_data |> select(-c(1:8)))

skew <- data.frame(Skewness=NA, p=NA, variable=NA)

for(var in variables) {
  tmp <- test.skew(sd_data$var) |>
    mutate(variable=var)
  
  skew <- bind_rows(skew, tmp)
}
str(sd_data)



test.skew(sd_data$CHLA)
