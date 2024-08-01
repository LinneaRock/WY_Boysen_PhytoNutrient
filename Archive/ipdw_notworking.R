
source('Data/CALL_DATA_LIB.R')

library(ipdw)
library(sf)
library(tidyverse)

pols <- st_read('C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')
st_crs(pols)
pols <- st_transform(pols, crs = 4326)

pnts <- read.csv('Data/RawData_WYDEQ/ChemPhysData_2002-2021.csv', 
                 fileEncoding="latin1") |> # doesn't matter which csv we call
  mutate(CollDate = as.Date(CollDate, format='%m/%d/%Y'),
         Year = year(CollDate)) |>
  filter(Year>=2020) |>
  dplyr::select(WaterbodyName, Latitude, Longitude) |>
  filter(grepl('Boysen', WaterbodyName)) |> # keep only Boysen 
  mutate(WaterbodyName = sub("^[^,]*,", "", WaterbodyName)) |> # shorten names since we know its Boysen
  distinct() |>
  mutate(WaterbodyName=trimws(WaterbodyName),
         lat=Latitude,
         lon=Longitude) |>
  mutate(WaterbodyName=factor(WaterbodyName, levels=c('Lacustrine Pelagic: Dam', 'East Shore','Cottonwood Creek Bay','Tough Creek Campground','Transitional Pelagic: Sand Mesa','Riverine Pelagic: Freemont 1','Fremont Bay'))) |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326) 

plot(pols)
plot(pnts)

costras <- costrasterGen(pnts, pols, extent = "pols",
                         projstr = projection(pols),resolution=0.0001)
plot(costras)
# insert contiguous barrier
costras[160:170, 1:80] <- 10000

st_crs(pnts) == st_crs(pols) 
st_crs(pols)==st_crs(costras)

# find average nearest neighbor
library(spatstat)

W              <- owin(range(c(st_bbox(pnts)["xmin"], st_bbox(pnts)["xmax"])),
                       range(c(st_bbox(pnts)["ymin"], st_bbox(pnts)["ymax"])))
kat.pp         <- ppp(st_coordinates(pnts)[,1], st_coordinates(pnts)[,2], window = W)
mean.neighdist <- mean(nndist(kat.pp))

# grid building
gridsize       <- mean.neighdist * 2
grainscale.fac <- gridsize / res(costras)[1]
gridras        <- aggregate(costras, fact = grainscale.fac)
gridpol        <- rasterToPolygons(gridras)
gridpol$value  <- row.names(gridpol)

# spatial join
fulldataset.over <- sf::st_join(pnts, st_as_sf(gridpol))

# grid selection
library(gdata)
set.seed(2)
gridlev <- unique(fulldataset.over$value)
for (i in seq_along(gridlev)) {
  activesub <- subset(fulldataset.over, fulldataset.over$value == gridlev[i])
  selectnum <- gdata::resample(seq_len(nrow(activesub)), 1)
  if (i == 1) {
    training <- activesub[selectnum, ]
  } else {
    training <- rbind(training, activesub[selectnum, ])
  }
}

validate             <- fulldataset.over[!(row.names(fulldataset.over) %in%
                                             row.names(training)), ]


plot(costras)
plot(st_geometry(training), add = TRUE)
plot(st_geometry(validate), col = "red", add = TRUE)

points <- rbind(training,validate)


raw_dat <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  bind_rows(BoysenPhyto_cat |> pivot_longer(cols=c(Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate), names_to = 'ShortName_Revised', values_to = 'ChemValue') |>
              left_join(BoysenChem |> dplyr::select(WaterbodyName, Latitude, Longitude) |> distinct())) |>
  dplyr::select(WaterbodyName, Year, month, Latitude, Longitude, ShortName_Revised, ChemValue) |>
  st_as_sf(coords=c('Longitude', 'Latitude'), crs=4326) |>
  st_join(points) 

plot(costras)
plot(st_geometry(raw_dat), add = TRUE)

raw_dat1 <- raw_dat |>
  filter(Year=='2023',
         month=='Jun',
         ShortName_Revised=='TN')

paramlist <- c('ChemValue')
final.ipdw <- ipdw(raw_dat1, costras, range = mean.neighdist * 10, paramlist,
                   overlapped = TRUE)
plot(final.ipdw)
?ipdw()
