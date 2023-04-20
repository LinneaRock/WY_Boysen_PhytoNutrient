#----------------------------------#
# Script to call libraries and data 
#----------------------------------#

library(tidyverse)
library(lubridate)


ChemPhys <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                     fileEncoding="latin1") |>
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate))
