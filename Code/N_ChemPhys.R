#-----------------------------------------------#
# Script to look at some basics of ChemPhys data
#-----------------------------------------------#

# 1. N ####
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
  distinct()

total_belowDet <- n_chemphys |>  
  mutate(belowDET_perc = BelowDetectionSite/TotalCollSite * 100) |>
  select(StationID, WaterbodyName, ShortName_Revised, belowDET_perc) |>
  distinct()
