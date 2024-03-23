# Example of plotting heatmaps of lake profile data ####

# 1. Load libraries and data ####
library(tidyverse)

CWC_profile <- BoysenProfile |>
  filter(WaterbodyName=='Cottonwood Creek Bay') |>
  rename(sampledate=CollDate,
         depth=depth_m)


observationDF=CWC_profile
# 2. Write function to interpolate between sampled points ####
## Simple vertical linear interpolation of water column concentrations function ####
# see more at https://github.com/hdugan/NTLlakeloads/tree/master
interpData <- function(observationDF, date, maxdepth) {
  a = observationDF %>% filter(sampledate == sampledate)
  if (sum(!is.na(a$param)) == 0) {
    print('nothing')
    return(NULL)
  }
  
  b = a %>% filter(!is.na(param))
  if (max(b$depth) < (maxdepth/2)) {
    print('too shallow')
    return(NULL)
  }
  
  yout = approx(x = a$depth, y = a$param, xout = xout, rule = 2)
  return(yout$y)
}

# 3. Run function on the data ####
# parameters
maxdepth<-max(CWC_profile$depth) # maxdepth of lowest sample, not necessarily depth of lake
param<-'temp_C'
xout=c(0,5,10,round(maxdepth)) 

usedates <- CWC_profile |>
  distinct(sampledate)


# run function
f <- lapply(X = usedates$sampledate, FUN = interpData, observationDF = CWC_profile, maxdepth = maxdepth)


f = as.data.frame(do.call(cbind, f))
names(f) = usedates$sampledate

# Bind list into dataframe
f2 = bind_cols(depth = c(0, 5, 10, 15, 20, 23.5),f) |>
  pivot_longer(-1, names_to = 'sampledate', values_to = 'var') |>
  arrange(sampledate,depth) |>
  mutate(sampledate = as.Date(sampledate))

# 4. Make a heatmap ####
ggplot() +
  guides(fill = guide_colorsteps(barheight = unit(3, "in"))) +
  geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
  geom_point(ME_profiles, mapping = aes(sampledate, depth), color = "#C4CFD0") +
  #theme_minimal() +
  scale_y_reverse() +
  scale_fill_viridis_d("Chloride
Concentration"~(mg~L^-1), option = "inferno", direction = -1) +
  labs(x = "", y = "Depth (m)") +
  theme(plot.caption = element_text(size = 10, hjust = 0), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  theme_bw()




ggplot(CWC_profile, aes(sampledate, depth, color = temp_C)) +
  geom_point() +
  scale_fill_viridis_c()


ggplot() +
  guides(fill = guide_colorsteps(barheight = unit(3, "in"))) +
  geom_contour_filled(CWC_profile, mapping = aes(x = sampledate, y = depth, z = temp_C), binwidth = 1) +
  geom_point(CWC_profile, mapping = aes(sampledate, depth), color = "#C4CFD0") +
  #theme_minimal() +
  scale_y_reverse() +
  scale_fill_viridis_d("Chloride
Concentration"~(mg~L^-1), option = "inferno", direction = -1) +
  labs(x = "", y = "Depth (m)") +
  theme(plot.caption = element_text(size = 10, hjust = 0), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  theme_bw()