
library(tidyverse)
library(dataRetrieval)


muddycrkparams <- whatNWISdata(siteNumber = "06258000")
fivemilecrkparams <- whatNWISdata(siteNumber = "06253000")
windrvrparams <- whatNWISdata(siteNumber = "06236100")
windrvroutletparams <- whatNWISdata(siteNumber = '06259000')
