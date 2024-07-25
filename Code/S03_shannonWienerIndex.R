#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Shannon-Wiener Index Figure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Shannon-Weiner Diversity Index; <1.5 low diversity, >2.5 high diversity'
# Calculation occurs in Data/CALL_DATA_LIB.R

source('Data/CALL_DATA_LIB.R')

# determine how different sampling locations are in terms of stability 
h <- aov(H~month, BoysenPhyto_A|> select(month, WaterbodyName, H) |> distinct())
tukey <- TukeyHSD(h)
library(multcompView)
cld <- multcompLetters4(h, tukey)
cld2 <- data.frame(letters = cld$month$Letters)
cld2$month <- rownames(cld2)
sig.letters <- cld2

# significance for plotting
sig.letters <- sig.letters |>
drop_na(letters) 

means <- left_join(BoysenPhyto_A |> 
                     select(month, WaterbodyName, H) |> 
                     distinct(), sig.letters) |>
  group_by(month, letters) |>
  summarise(max.result = max(H, na.rm = TRUE)) |>
  distinct()

# make plot
BoysenPhyto_A |>
  select(month, WaterbodyName, H) |> 
  distinct() |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  ggplot() +
  geom_boxplot(aes(month, H)) +
  geom_jitter(aes(month, H, color=WaterbodyName),alpha=0.5) +
  scale_color_viridis_d('', option='turbo') +
  geom_text(means, mapping=aes(month, 
                               max.result+0.5, label = letters), 
            size=4) +
  labs(x='',y='Shannon-Weiner diversity index') +
  theme_minimal() 
ggsave('Figures/H_diversity.png',height=4.5,width=6.5,units='in',dpi=1200)
