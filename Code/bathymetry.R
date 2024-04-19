library(terra)
library(tidyterra)



#If using the tif DEM
boysen_bathy <- rast("C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/BoysenRaster.tif")

#To plot contours easily all you have to do is:
contour(boysen_bathy) #this only works with the raster not the vector

#If using the contour lines shapefile
boysen_bathy <- vect("C:/Users/linne/OneDrive - University of Wyoming/Data/Spatial_Data/Boysen/Boysen_Bathy/Boysen_Contours/Boysen_Contours/All_Boysen_Contours.shp")
#You can plot either of these easily using base plot or use the tidyterra package for ggplot


#As far as cross sectional analysis you'll probably want to use the raster from the DEM tif. Terra has awesome documentation that gives a great idea about the tools available to you, but if you want to talk more in person about what you want to do I'd be happy to.