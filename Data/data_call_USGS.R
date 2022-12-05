
library(tidyverse)
library(dataRetrieval)

# What data are available at each site? ####
muddycrkparams <- whatNWISdata(siteNumber = "06258000")
fivemilecrkparams <- whatNWISdata(siteNumber = "06253000")
windrvrparams <- whatNWISdata(siteNumber = "06236100")
windrvroutletparams <- whatNWISdata(siteNumber = '06259000')

# Load site data ####
sitedata <- readNWISsite(siteNumbers = c("06258000", "06253000", "06236100", '06259000'))
                         
# Load parameters we want ####
alldat <- readNWISqw(siteNumbers = c("06258000", "06253000", "06236100", '06259000'),
                          parameterCd = c('00010', '00060', '000095', '00400', '00600', '00605', '00608', '00631', '00665','00671'),
                          startDate = '2002-01-01') #|>
  # mutate(parameter = case_when(param_cd == '00600' ~ 'TN_mgL', 
  #                         param_cd == '00605' ~ 'TON_mgL', 
  #                         param_cd == '00613' ~ 'Wind River Inlet', 
  #                         param_cd == '00631' ~ 'Wind River Outlet',
  #                         param_cd == '00665',
  #                         param_cd == '00671'))
