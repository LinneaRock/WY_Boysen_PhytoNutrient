#-----------------------------------------------#
# Exploratory plots for April 2023 meeting
#-----------------------------------------------#

source('Data/CALL_DATA_LIB.R')

params <- unique(ChemPhys$ShortName_Revised) # which parameters are available in this dataset? 
params

# 1. Prep data for plotting ####
BoysenChemPhys <- ChemPhys |>
  select(-ChemSampID) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  mutate(ChemUnits = sub('/l', '/L', ChemUnits)) |>
  # unite(col=Param, ShortName_Revised, ChemUnits, sep='_') |> # combine name and units
  # mutate(Param= gsub('[^[:alnum:]]+', '_', Param)) |> # fix for column header (for now)
  select(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, BelowDet, ChemValue, ChemUnits, SampleDepth) |> # keep just what is interesting now
  filter(ShortName_Revised %in% c("Phosphorus as P (total)","Total Nitrogen (unfiltered)","Nitrate as N", "Orthophosphate as P (total)","Nitrate plus Nitrite as N","Total Ammonia as N","Chlorophyll a (phytoplankton)","Total Kjeldahl Nitrogen (unfiltered)","Nitrite as N")) |> # for now just keep at the nutrients and chlorophyll
  distinct() |>  # remove duplicates 
  group_by(StationID, WaterbodyName, Latitude, Longitude, CollDate, ShortName_Revised, ChemUnits, SampleDepth) |>
  mutate(ChemValue = mean(ChemValue)) |> # take average of replicates
  ungroup() |>
  distinct()

## how many of each parameter are below detection at each location? ####
n_belowDet <- BoysenChemPhys |>
  group_by(WaterbodyName, ShortName_Revised) |>
  add_count() |>
  mutate(range = paste(min(ChemValue), max(ChemValue), sep='-')) |>
  ungroup() |>
  select(WaterbodyName, ShortName_Revised, BelowDet, n, range) |>
  group_by(WaterbodyName, ShortName_Revised, n, range) |>
  summarise(belowDet = sum(BelowDet),
            percLow = belowDet/n *100) |>
  distinct()

# 2. Basic timeseries! ####
params <- unique(BoysenChemPhys$ShortName_Revised)

for(p in 1:length(params)) {

    ggplot(BoysenChemPhys |> filter(ShortName_Revised==params[p])) + 
      geom_point(aes(CollDate, SampleDepth, color=ChemValue)) +
      labs(y='Sample depth (m)', x='', title=params[p]) +
      scale_y_reverse() +
      facet_wrap(~WaterbodyName, scales='free_y') +
      scale_color_viridis_c(paste(params[p],unique((BoysenChemPhys |>
                                                      filter(ShortName_Revised==params[p]))$ChemUnits), 
                                  sep = ' ')) +
      theme_classic() 
    ggsave(paste0('Figures/Boysen_explore/',params[p],'_depth.png'),width=10, height=8, units='in')
  
    ggplot(BoysenChemPhys |> filter(ShortName_Revised==params[p],
                                    SampleDepth <= 3)) + 
      geom_point(aes(CollDate, ChemValue)) +
      labs(y=paste(params[p],unique((BoysenChemPhys |> 
                                       filter(ShortName_Revised==params[p]))$ChemUnits), 
                   sep = ' '), x='') +
      scale_y_reverse() +
      facet_wrap(~WaterbodyName, scales='free_y') +
      theme_classic() 
    ggsave(paste0('Figures/Boysen_explore/',params[p],'_epi.png'),width=10, height=8, units='in')
    
}




