## heatmaps of lake profile data ####
library(patchwork)
# some of these are not working for certain locations/seasons


source('Data/CALL_DATA_LIB.R')

BoysenProfile$WaterbodyName <- gsub(":.*", "", BoysenProfile$WaterbodyName)


WaterbodyName <- unique(BoysenProfile$WaterbodyName)



# WATER TEMPERATURE ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring

interp_plot_temp <- function(location) {

 ## 1.. Filter data ####
  # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
boy_pro <- BoysenProfile |>
  rename(sampledate=CollDate,
         depth=depth_m,
         param=temp_C) |>
  filter(WaterbodyName==location) |>
  select(WaterbodyName, sampledate, depth, param) |>
  mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', NA)) |>
  mutate(summer=ifelse(year(sampledate)==2021, 'Summer 2021',
                              ifelse(year(sampledate)==2022, 'Summer 2022',
                                     ifelse(year(sampledate)==2023,'Summer 2023',summer)))) 

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
  # get rid of this becuase we have some shallow sites - did checks and plots look fine with/without the depth justification
  # if (max(b$depth) < (maxdepth/2)) {
  #   print('too shallow')
  #   return(NULL)
  # }
  
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
  mutate(sampledate = as.Date(sampledate)) |>
  left_join(boy_pro |> select(sampledate,summer))

# 4. Make a heatmap ####
ggplot() +
  guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
  geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
  geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
  scale_fill_viridis_d("Temperature"~(degree*C), option='inferno') +
  scale_y_reverse() +
  labs(x = "", y = "Depth (m)") +
  theme(plot.caption = element_text(size = 10, hjust = 0), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  theme_minimal() +
  facet_wrap(~summer, ncol = 2, scale='free_x')

}

patched_plot <- interp_plot_temp(location=name)
patched_plot + plot_annotation(title = name) 
ggsave(paste0('Figures/Profiles/Temp/',name,'.png'),height=4.5,width=6.5,units='in',dpi=1200)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch

}







# pH ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_ph <- function(location) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=pH) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', NA)) |>
        mutate(summer=ifelse(year(sampledate)==2021, 'Summer 2021',
                             ifelse(year(sampledate)==2022, 'Summer 2022',
                                    ifelse(year(sampledate)==2023,'Summer 2023',summer)))) 
      
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
        # get rid of this becuase we have some shallow sites - did checks and plots look fine with/without the depth justification
        # if (max(b$depth) < (maxdepth/2)) {
        #   print('too shallow')
        #   return(NULL)
        # }
        
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
        mutate(sampledate = as.Date(sampledate)) |>
        left_join(boy_pro |> select(sampledate,summer))
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d('pH', option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)") +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal() +
        facet_wrap(~summer, ncol = 2, scale='free_x')
      
    }
    
    patched_plot <- interp_plot_ph(location=name)
    patched_plot + plot_annotation(title = name) 
    ggsave(paste0('Figures/Profiles/pH/',name,'.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}










# SpC ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_spc <- function(location) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=cond_uScm) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', NA)) |>
        mutate(summer=ifelse(year(sampledate)==2021, 'Summer 2021',
                             ifelse(year(sampledate)==2022, 'Summer 2022',
                                    ifelse(year(sampledate)==2023,'Summer 2023',summer)))) 
      
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
        # get rid of this becuase we have some shallow sites - did checks and plots look fine with/without the depth justification
        # if (max(b$depth) < (maxdepth/2)) {
        #   print('too shallow')
        #   return(NULL)
        # }
        
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
        mutate(sampledate = as.Date(sampledate))|>
        left_join(boy_pro |> select(sampledate,summer) |> distinct())
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var)) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d('SpC '~(µS~cm^-1), option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)") +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal() +
        facet_wrap(~summer, ncol = 2, scale='free_x')
      
    }
    
    patched_plot <- interp_plot_spc(location=name)
    patched_plot + plot_annotation(title = name) 
    ggsave(paste0('Figures/Profiles/SpC/',name,'.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}








# DO ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_do <- function(location) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=DO_mgL) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', NA)) |>
        mutate(summer=ifelse(year(sampledate)==2021, 'Summer 2021',
                             ifelse(year(sampledate)==2022, 'Summer 2022',
                                    ifelse(year(sampledate)==2023,'Summer 2023',summer)))) 
      
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
        # get rid of this becuase we have some shallow sites - did checks and plots look fine with/without the depth justification
        # if (max(b$depth) < (maxdepth/2)) {
        #   print('too shallow')
        #   return(NULL)
        # }
        
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
        mutate(sampledate = as.Date(sampledate))|>
        left_join(boy_pro |> select(sampledate,summer) |> distinct())
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d('DO'~(mg~L^-1), option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)") +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal() +
        facet_wrap(~summer, ncol = 2, scale='free_x')
      
    }
    
    patched_plot <- interp_plot_do(location=name)
    patched_plot + plot_annotation(title = name) 
    ggsave(paste0('Figures/Profiles/DO/',name,'.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}


# supsected faulty readings in LPD on 2021-09-15
lpd_rm <- BoysenProfile |>
  filter(WaterbodyName=='Lacustrine Pelagic: Dam') |>
  mutate(DO_mgL = ifelse(CollDate==as.Date('2021-09-15'), NA, DO_mgL)) |>
  rename(sampledate=CollDate,
         depth=depth_m,
         param=DO_mgL) |>
  select(WaterbodyName, sampledate, depth, param) |>
  mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', NA)) |>
  mutate(summer=ifelse(year(sampledate)==2021, 'Summer 2021',
                       ifelse(year(sampledate)==2022, 'Summer 2022',
                              ifelse(year(sampledate)==2023,'Summer 2023',summer)))) |>
  drop_na()

# run interpdata for function then:

# parameters
maxdepth <- max(lpd_rm $depth) # maxdepth of lowest sample, not necessarily depth of lake
usedates <- lpd_rm  |>
  distinct(sampledate)
xout<-seq(0,maxdepth+1,by=1)


# run function
f <- lapply(X = usedates$sampledate, FUN = interpData, observationDF = lpd_rm, maxdepth = maxdepth)


f = as.data.frame(do.call(cbind, f))
names(f) = usedates$sampledate

# Bind list into dataframe
f2 = bind_cols(depth = xout,f) |>
  pivot_longer(-1, names_to = 'sampledate', values_to = 'var') |>
  arrange(sampledate,depth) |>
  mutate(sampledate = as.Date(sampledate))|>
  left_join(lpd_rm  |> select(sampledate,summer) |> distinct())

# 4. Make a heatmap ####
ggplot() +
  guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
  geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var), binwidth = 1) +
  geom_point(lpd_rm, mapping = aes(sampledate, depth), color = "#C4CFD0") +
  scale_fill_viridis_d('DO'~(mg~L^-1), option='inferno') +
  scale_y_reverse() +
  labs(x = "", y = "Depth (m)", title='Lacustrine Pelagic - removed faulty reading') +
  theme(plot.caption = element_text(size = 10, hjust = 0), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  theme_minimal() +
  facet_wrap(~summer, ncol = 2, scale='free_x')


ggsave(paste0('Figures/Profiles/DO/LPD_remove09152021.png'),height=4.5,width=6.5,units='in',dpi=1200)










# ORP ####

for(name in WaterbodyName) {
  tryCatch({ # continues running loop, but gives error message where problem(s) are occurring
    
    interp_plot_orp <- function(location) {
      
      ## 1.. Filter data ####
      # can't figure out how to add column in the function so we'll have separate functions for each parameter :( 
      boy_pro <- BoysenProfile |>
        rename(sampledate=CollDate,
               depth=depth_m,
               param=ORP) |>
        filter(WaterbodyName==location) |>
        select(WaterbodyName, sampledate, depth, param) |>
        mutate(summer=ifelse(year(sampledate)==2020, 'Summer 2020', NA)) |>
        mutate(summer=ifelse(year(sampledate)==2021, 'Summer 2021',
                             ifelse(year(sampledate)==2022, 'Summer 2022',
                                    ifelse(year(sampledate)==2023,'Summer 2023',summer)))) 
      
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
        # get rid of this becuase we have some shallow sites - did checks and plots look fine with/without the depth justification
        # if (max(b$depth) < (maxdepth/2)) {
        #   print('too shallow')
        #   return(NULL)
        # }
        
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
        mutate(sampledate = as.Date(sampledate)) |>
        left_join(boy_pro |> select(sampledate,summer))
      
      # 4. Make a heatmap ####
      ggplot() +
        guides(fill = guide_colorsteps(barheight = unit(2, "in"))) +
        geom_contour_filled(f2, mapping = aes(x = sampledate, y = depth, z = var)) +
        geom_point(boy_pro, mapping = aes(sampledate, depth), color = "#C4CFD0") +
        scale_fill_viridis_d('ORP', option='inferno') +
        scale_y_reverse() +
        labs(x = "", y = "Depth (m)") +
        theme(plot.caption = element_text(size = 10, hjust = 0), 
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 10)) +
        theme_minimal() +
        facet_wrap(~summer, ncol = 2, scale='free_x')
      
    }
    
    patched_plot <- interp_plot_orp(location=name)
    patched_plot + plot_annotation(title = name) 
    ggsave(paste0('Figures/Profiles/ORP/',name,'.png'),height=4.5,width=6.5,units='in',dpi=1200)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # end try catch
  
}







