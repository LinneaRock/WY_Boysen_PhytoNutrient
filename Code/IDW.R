#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# inverse distance weighted interpolation across Boysen surface
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')
library(gstats)
library(sp)
library(sf)


# 2. read in lake shapefile ####
lake_shapefile <- st_read('C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen Shapefile/Boysen_Shape.shp')
st_crs(lake_shapefile)
lake_shapefile <- st_transform(lake_shapefile, crs = 3738) # https://epsg.io/3738

ggplot()+
  geom_sf(lake_shapefile, mapping=aes())

# 3. create grid ####
#make a grid across shp
#This grid can be saved and read in to streamline future grid usage
#cellsize are based on CRS, make sure to use a UTM system
grid <- lake_shapefile  %>%
  st_make_grid(cellsize = 1000, what = "centers", square = T) %>% # grid of points
  st_intersection(lake_shapefile )


#They look okay?
ggplot(grid)+
  geom_sf(size = 0.1)


# 4. Format data ####
raw_dat <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  bind_rows(BoysenPhyto_cat |> pivot_longer(cols=c(Diatom, `Green algae`,  Cyanobacteria, Dinoflagellate, `Golden algae`, Flagellate), names_to = 'ShortName_Revised', values_to = 'ChemValue') |>
              left_join(BoysenChem |> dplyr::select(WaterbodyName, Latitude, Longitude) |> distinct())) |>
  dplyr::select(WaterbodyName, Year, month, Latitude, Longitude, ShortName_Revised, ChemValue) |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326)

raw_dat <- st_transform(raw_dat, 3738)

#make sure data is there
ggplot()+
  geom_sf(data = lake_shapefile)+
  geom_sf(data = grid)+
  geom_sf(data = raw_dat, aes(), color='red') 


# 5. IDW ####

# run over each sampling date and for each parameter
months <- c('May','Jun','Jul','Aug','Sep','Oct')
years <- c('2020','2021','2022','2023')
vars <- as.vector(as.data.frame(raw_dat) |> dplyr::select(ShortName_Revised) |> distinct())[["ShortName_Revised"]]


for(y in years) {
  for(m in months) {
    for(v in vars) {
      tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
      
      tmp <- raw_dat |>
        filter(month==m,
               Year==y,
               ShortName_Revised==v) 
      
      idw_dat <-gstat::idw(formula = tmp$ChemValue ~ 1,
                           locations = tmp,
                           newdata = grid,
                           idp = 3)
      
      #plot the idw data
    p <- ggplot()+
        ggtitle(paste0(y,': ', m, '-', v))+
        geom_sf(data=lake_shapefile, aes()) +
        geom_sf(data = idw_dat, aes(color = var1.pred))+
        theme_bw() +
        theme(axis.text.x=element_text(angle=45))
      
      ggsave(filename=paste0('Figures/IDW/',v,'/',y,m,'.png'),plot=p,height=4.5,width=6.5,units='in',dpi=1200)
      
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
      
    }
  }
}











