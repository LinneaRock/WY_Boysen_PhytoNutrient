## heatmaps of lake profile data ####

# some of these are not working for certain locations/seasons


source('Data/CALL_DATA_LIB.R')

BoysenProfile$WaterbodyName <- gsub(":.*", "", BoysenProfile$WaterbodyName)


WaterbodyName <- unique(BoysenProfile$WaterbodyName)



# WATER TEMPERATURE ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring

interp_plot_temp <- function(location,season, legend) {

 ## 1.. Filter data ####
  # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
boy_pro <- BoysenProfile |>
  rename(sampledate=CollDate,
         depth=depth_m,
         param=temp_C) |>
  filter(WaterbodyName==location) |>
  select(WaterbodyName, sampledate, depth, param) |>
  mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', 'Summer 2021')) |>
  filter(summer==season)

# 2. Write function to interpolate between sampled points ####
## Simple vertical linear interpolation of water column concentrations function ####
# see more at https://github.com/hdugan/NTLlakeloads/tree/master
interpData <- function(observationDF, date, maxdepth) {
  a = observationDF %>% filter(sampledate == date)
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
maxdepth <- max(boy_pro$depth) # maxdepth of lowest sample, not necessarily depth of lake
usedates <- boy_pro |>
  distinct(sampledate)
xout<-seq(0,maxdepth+1,by=1)


# run function
f <- lapply(X = usedates$sampledate, FUN = interpData, observationDF = boy_pro, maxdepth = maxdepth)


f = as.data.frame(do.call(cbind, f))
names(f) = usedates$sampledate

# Bind list into dataframe
f2 = bind_cols(depth = xout,f) |>
  pivot_longer(-1, names_to = 'sampledate', values_to = 'var') |>
  arrange(sampledate,depth) |>
  mutate(sampledate = as.Date(sampledate))

# 4. Make a heatmap ####
ggplot() +
  guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
  geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
  geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
  scale_fill_viridis_d(legend, option='inferno') +
  scale_y_reverse() +
  labs(x = "", y = "Depth (m)", title=season) +
  theme(plot.caption = element_text(size = 10, hjust = 0), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  theme_minimal()

}



tmp1 <- interp_plot_temp(location=name,season='Summer 2020',legend=expression("Temperature"~(degree*C)))

tmp2 <- interp_plot_temp(location=name,season='Summer 2021',legend=expression("Temperature"~(degree*C)))



tmp1 +
  plot_annotation(title = name) 

ggsave(paste0('Figures/Profiles/temp/',name,'a.png'),height=4.5,width=6.5,units='in',dpi=1200)

tmp2 +
  plot_annotation(title = name) 

ggsave(paste0('Figures/Profiles/temp/',name,'b.png'),height=4.5,width=6.5,units='in',dpi=1200)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch

}







# pH ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_ph <- function(location,season, legend) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=pH) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', 'Summer 2021')) |>
        filter(summer==season)
      
      # 2. Write function to interpolate between sampled points ####
      ## Simple vertical linear interpolation of water column concentrations function ####
      # see more at https://github.com/hdugan/NTLlakeloads/tree/master
      interpData <- function(observationDF, date, maxdepth) {
        a = observationDF %>% filter(sampledate == date)
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
      maxdepth <- max(boy_pro$depth) # maxdepth of lowest sample, not necessarily depth of lake
      usedates <- boy_pro |>
        distinct(sampledate)
      xout<-seq(0,maxdepth+1,by=1)
      
      
      # run function
      f <- lapply(X = usedates$sampledate, FUN = interpData, observationDF = boy_pro, maxdepth = maxdepth)
      
      
      f = as.data.frame(do.call(cbind, f))
      names(f) = usedates$sampledate
      
      # Bind list into dataframe
      f2 = bind_cols(depth = xout,f) |>
        pivot_longer(-1, names_to = 'sampledate', values_to = 'var') |>
        arrange(sampledate,depth) |>
        mutate(sampledate = as.Date(sampledate))
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d(legend, option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)", title=season) +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal()
      
    }
    
    
    
    tmp1 <- interp_plot_ph(location=name,season='Summer 2020',legend='pH')
    
    tmp2 <- interp_plot_ph(location=name,season='Summer 2021',legend='pH')
    
    
    
    tmp1 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/ph/',name,'a.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
    tmp2 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/ph/',name,'b.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}











# SpC ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_cond <- function(location,season, legend) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=cond_ugL) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', 'Summer 2021')) |>
        filter(summer==season)
      
      # 2. Write function to interpolate between sampled points ####
      ## Simple vertical linear interpolation of water column concentrations function ####
      # see more at https://github.com/hdugan/NTLlakeloads/tree/master
      interpData <- function(observationDF, date, maxdepth) {
        a = observationDF %>% filter(sampledate == date)
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
      maxdepth <- max(boy_pro$depth) # maxdepth of lowest sample, not necessarily depth of lake
      usedates <- boy_pro |>
        distinct(sampledate)
      xout<-seq(0,maxdepth+1,by=1)
      
      
      # run function
      f <- lapply(X = usedates$sampledate, FUN = interpData, observationDF = boy_pro, maxdepth = maxdepth)
      
      
      f = as.data.frame(do.call(cbind, f))
      names(f) = usedates$sampledate
      
      # Bind list into dataframe
      f2 = bind_cols(depth = xout,f) |>
        pivot_longer(-1, names_to = 'sampledate', values_to = 'var') |>
        arrange(sampledate,depth) |>
        mutate(sampledate = as.Date(sampledate))
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d(legend, option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)", title=season) +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal()
      
    }
    
    
    
    tmp1 <- interp_plot_cond(location=name,season='Summer 2020',legend='SpC'~(µS~cm^-1))
    
    tmp2 <- interp_plot_cond(location=name,season='Summer 2021',legend='SpC '~(µS~cm^-1))
    
    
    
    tmp1 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/cond/',name,'a.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
    tmp2 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/cond/',name,'b.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}







# DO ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_do <- function(location,season, legend) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=DO_mgL) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', 'Summer 2021')) |>
        filter(summer==season)
      
      # 2. Write function to interpolate between sampled points ####
      ## Simple vertical linear interpolation of water column concentrations function ####
      # see more at https://github.com/hdugan/NTLlakeloads/tree/master
      interpData <- function(observationDF, date, maxdepth) {
        a = observationDF %>% filter(sampledate == date)
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
      maxdepth <- max(boy_pro$depth) # maxdepth of lowest sample, not necessarily depth of lake
      usedates <- boy_pro |>
        distinct(sampledate)
      xout<-seq(0,maxdepth+1,by=1)
      
      
      # run function
      f <- lapply(X = usedates$sampledate, FUN = interpData, observationDF = boy_pro, maxdepth = maxdepth)
      
      
      f = as.data.frame(do.call(cbind, f))
      names(f) = usedates$sampledate
      
      # Bind list into dataframe
      f2 = bind_cols(depth = xout,f) |>
        pivot_longer(-1, names_to = 'sampledate', values_to = 'var') |>
        arrange(sampledate,depth) |>
        mutate(sampledate = as.Date(sampledate))
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d(legend, option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)", title=season) +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal()
      
    }
    
    
    
    tmp1 <- interp_plot_do(location=name,season='Summer 2020',legend='DO'~(mg~L^-1))
    
    tmp2 <- interp_plot_do(location=name,season='Summer 2021',legend='DO'~(mg~L^-1))
    
    
    
    tmp1 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/DO/',name,'a.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
    tmp2 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/DO/',name,'b.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}
















# ORP ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_orp <- function(location,season, legend) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=ORP) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', 'Summer 2021')) |>
        filter(summer==season)
      
      # 2. Write function to interpolate between sampled points ####
      ## Simple vertical linear interpolation of water column concentrations function ####
      # see more at https://github.com/hdugan/NTLlakeloads/tree/master
      interpData <- function(observationDF, date, maxdepth) {
        a = observationDF %>% filter(sampledate == date)
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
      maxdepth <- max(boy_pro$depth) # maxdepth of lowest sample, not necessarily depth of lake
      usedates <- boy_pro |>
        distinct(sampledate)
      xout<-seq(0,maxdepth+1,by=1)
      
      
      # run function
      f <- lapply(X = usedates$sampledate, FUN = interpData, observationDF = boy_pro, maxdepth = maxdepth)
      
      
      f = as.data.frame(do.call(cbind, f))
      names(f) = usedates$sampledate
      
      # Bind list into dataframe
      f2 = bind_cols(depth = xout,f) |>
        pivot_longer(-1, names_to = 'sampledate', values_to = 'var') |>
        arrange(sampledate,depth) |>
        mutate(sampledate = as.Date(sampledate))
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d(legend, option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)", title=season) +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal()
      
    }
    
    
    
    tmp1 <- interp_plot_orp(location=name,season='Summer 2020',legend='ORP')
    
    tmp2 <- interp_plot_orp(location=name,season='Summer 2021',legend='ORP')
    
    
    
    tmp1 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/ORP/',name,'a.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
    tmp2 +
      plot_annotation(title = name) 
    
    ggsave(paste0('Figures/Profiles/ORP/',name,'b.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}








