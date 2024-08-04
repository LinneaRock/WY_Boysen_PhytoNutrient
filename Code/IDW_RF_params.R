#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# inverse distance weighted interpolation across Boysen surface
# top 4 random forest predictors 2022 only, including mapping
# cyanotoxin locations 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. load libraries, data ####
source('Data/CALL_DATA_LIB.R')
library(gstat)
library(sp)
library(sf)


# 2. read in lake shapefile ####
# shapefile too big to be stored on github
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
# keep top 4 predictors and 2022 data only for example
raw_dat <- BoysenNutrient |>
  bind_rows(BoysenChem) |>
  dplyr::select(WaterbodyName, Year, month, Latitude, Longitude, ShortName_Revised, ChemValue) |>
  filter(ShortName_Revised %in% c('Secchi', 'SpC', 'DO', 'TP'),
         Year==2022) |>
  st_as_sf(coords=c('Longitude','Latitude'), crs=4326)

raw_dat <- st_transform(raw_dat, 3738)


# format cyanotoxin data for plotting locations
cyano_prep <- cyanotoxin |>
  mutate(month = month(CollDate,label=TRUE,abbr=TRUE),
         Year=year(CollDate)) |>
  dplyr::select(-CollDate, -Advisory) |>
  mutate(toxinpresent=as.character(toxinpresent)) |>
  mutate(toxinpresent=ifelse(toxinpresent=='N',0,1)) |>
  filter(toxinpresent != 0) |>
  distinct() |>
  filter(Year == 2022) 

# Create a new row with 'May' in the month column and NA for other columns
new_row <- data.frame(matrix(ncol = ncol(cyano_prep), nrow = 1))
colnames(new_row) <- colnames(cyano_prep)
new_row$month <- 'Jun'


# Bind the new row to the existing data frame
cyano_prep <- rbind(cyano_prep, new_row) 


# Create a new row with 'May' in the month column and NA for other columns
new_row <- data.frame(matrix(ncol = ncol(cyano_prep), nrow = 1))
colnames(new_row) <- colnames(cyano_prep)
new_row$month <- 'May'

# Bind the new row to the existing data frame
cyano_prep <- rbind(cyano_prep, new_row) |>
  st_as_sf(coords=c('Long','Lat'), crs=4326)

cyano_dat <- st_transform(cyano_prep, 3738)

#make sure data is there
ggplot()+
  geom_sf(data = lake_shapefile)+
  geom_sf(data = grid)+
  geom_sf(data = raw_dat, aes(), color='red') +
  geom_sf(data = cyano_dat, aes(), color='#999933')


# 5. IDW ####

# run over each sampling date and for each parameter
months <- c('May','Jun','Jul','Aug','Sep','Oct')
vars <- as.vector(as.data.frame(raw_dat) |> dplyr::select(ShortName_Revised) |> distinct())[["ShortName_Revised"]]

idw_dat <- data.frame(var1.pred=NA,
                         var1.var=NA,
                      month=NA) 

idw_dat$geometry <- st_sfc(NA)
str(idw_dat)
idw_dat <- st_as_sf(idw_dat, crs=3738)

for(v in vars) {
  
      tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
        
        tmpa <- raw_dat |>
          filter(ShortName_Revised==v)
        
        for(m in months) {
        tmp <- tmpa |>
          filter(month==m ) 
        
        idw_dat_tmp <-gstat::idw(formula = tmp$ChemValue ~ 1,
                             locations = tmp,
                             newdata = grid,
                             idp = 3)
      }
       
        
         if (st_crs(idw_dat_tmp) != st_crs(idw_dat)) {
          idw_dat <- st_transform(idw_dat, st_crs(idw_dat_tmp))
        }
        
        
        idw_dat <- idw_dat |>
          bind_rows(idw_dat_tmp) |>
          mutate(month=m) |>
          drop_na(month) |>
          st_as_sf()
        
        # get minimum and maximum values for legend
        minlim <- (raw_dat |>
                     filter(ShortName_Revised == v) |>
                     summarise(min=min(ChemValue)))[[1]]
        
        maxlim <- (raw_dat |>
                     filter(ShortName_Revised == v) |>
                     summarise(max=max(ChemValue)))[[1]]
        
        # get cyano locations for that month
        # plot_loc <- cyano_dat |>
        #   filter(month==m)
        
        
        #plot the idw data
        p <- ggplot()+
          #ggtitle(paste0('2022: ', m))+
          geom_sf(data=lake_shapefile, aes()) +
          geom_sf(data = idw_dat, aes(color = var1.pred))+
          geom_sf(data=cyano_dat, aes(), color='#999933', size=2.5) +
          scale_color_gradient(paste0(v), high='grey98', low='grey2',limits=c(minlim, maxlim)) +
          theme_bw() +
          theme(axis.text.x=element_text(angle=45),
                axis.text.y=element_text(angle=45)) +
          facet_wrap(~month, nrow=1)
        
        ggsave(filename=paste0('Figures/IDW/2022example/',v,'/facetedToxin.png'),plot=p,height=4.5,width=6.5,units='in',dpi=1200)
        
        #above plots conservatively and won't include all months because they don't exist in the cyano dataset. So make another one without cyano data and manually fix plots in powerpoint :/
        p2 <- ggplot()+
         # ggtitle(paste0('2022: ', m))+
          geom_sf(data=lake_shapefile, aes()) +
          geom_sf(data = idw_dat, aes(color = var1.pred))+
          scale_color_gradient(paste0(v), high='grey98', low='grey2',limits=c(minlim, maxlim)) +
          theme_bw() +
          theme(axis.text.x=element_text(angle=45),
                axis.text.y=element_text(angle=45)) +
          facet_wrap(~month, nrow=1)
        
        ggsave(filename=paste0('Figures/IDW/2022example/',v,'/NOtoxins.png'),plot=p2,height=4.5,width=6.5,units='in',dpi=1200)
        
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
      
   
  }



