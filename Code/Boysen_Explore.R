#-----------------------------------------------#
# Exploratory plots for April 2023 meeting
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')


# 1. Prep data for plotting ####
BoysenChemPhys <- ChemPhys |>
  select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  #select(1:2,4:7, 12:16, 20, 23) |> # keep just what is interesting now
  distinct() |> # remove duplicates 
  group_by(StationID, WaterbodyName, CollDate, Param, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct() |>
  # pivot_wider( names_from = 'Param', values_from = 'ChemValue') |> # make each param a column for plotting
  # rename(# most names need to be fixed
  #   Chla_Phyto_ugL = Chlorophyll_a_phytoplankton_µg_l,
  #   Cond = Conductance_µS_cm,
  #   NO3_NO2_N_mgL = Nitrate_plus_Nitrite_as_N_mg_l,
  #   DO_mgL = DO_mg_L_mg_l,
  #   pH = pH_SU,
  #   NH4_N_mgL = Total_Ammonia_as_N_mg_l,
  #   TP_mgL = Phosphorus_as_P_total_mg_l,
  #   TN_mgL = Total_Nitrogen_unfiltered_mg_l,
  #   Ortho_P_mgL = Orthophosphate_as_P_total_mg_l) |>
  filter(BelowDet == 0) # remove low values for now
  

BoysenEpiChemPhys <- BoysenChemPhys |>
  filter(SampleDepth <= 3 | SampleDepth <= Secchi_Depth_m)

ggplot(BoysenEpiChemPhys)
  
