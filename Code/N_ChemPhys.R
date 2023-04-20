#-----------------------------------------------#
# Script to look at some basics of ChemPhys data
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')

# 1. N below detection ####
n_chemphys <- ChemPhys |>
  select(StationID, WaterbodyName, ShortName_Revised, Year, BelowDet) |>
  group_by(StationID, WaterbodyName, ShortName_Revised) |>
  add_count() |>
  rename(TotalCollSite = n) |>
  mutate(BelowDetectionSite = sum(BelowDet)) |>
  ungroup() |>
  group_by(StationID, WaterbodyName, ShortName_Revised, Year) |>
  add_count() |>
  rename(TotalCollYear = n) |>
  mutate(BelowDetectionYear = sum(BelowDet)) |>
  distinct()|>
  ungroup()

total_belowDet <- n_chemphys |>  
  mutate(belowDET_perc = BelowDetectionSite/TotalCollSite * 100) |>
  select(StationID, WaterbodyName, ShortName_Revised, belowDET_perc) |>
  distinct()
