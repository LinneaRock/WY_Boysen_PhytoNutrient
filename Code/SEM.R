#--------------------------------------------------------------------------#
# Structural equation modeling to predict high percentages of cyanobacteria
#--------------------------------------------------------------------------#

source('Data/CALL_DATA_LIB.R')


psych::pairs.panels(sd_data |> select(-c(1:8)))
