#----------------------------------#
# Script to call libraries and data 
#----------------------------------#

library(tidyverse)
library(lubridate)


ChemPhys <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                     fileEncoding="latin1") |>
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate))

BoysenTribs <- read.csv('Data/BoysenTribs.csv') |>
  mutate(CollDate = as.Date(CollDate, format='%Y-%m-%d'),
         Year = year(CollDate))


library(readr)
Phytoplankton_2013_2021 <- read_csv("Data/RawData_WYDEQ/Phytoplankton_2013-2021.csv", 
                                    skip = 5) |>
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate))

