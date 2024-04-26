#------------------------------------------------#
# Script to look at conditions when toxin present 
#------------------------------------------------#

source('Data/CALL_DATA_LIB.R')

nutrient_summarise <- BoysenNutrient |>
  group_by(CollDate, Year, ShortName_Revised) |>
  summarise(mean = mean(ChemValue),
            max = max(ChemValue),
            min = min(ChemValue))

chem_summarise <- BoysenChem |>
  group_by(CollDate, Year, ShortName_Revised) |>
  summarise(mean = mean(ChemValue),
            max = max(ChemValue),
            min = min(ChemValue))

data_summarise <- bind_rows(nutrient_summarise, chem_summarise)


ggplot(data_summarise) +
  geom_point(aes(CollDate, mean)) +
  geom_errorbar(aes(CollDate, mean, ymin=min, ymax=max, width=0.2)) +
  geom_point(cyanotoxin |> filter(toxinpresent=='Y'), mapping=aes(CollDate, 0), shape=3, color='red4') +
  theme_bw() +
  labs(x='',y='Reservoir average and range') +
  facet_wrap(~ShortName_Revised, scales='free',nrow=5) 


