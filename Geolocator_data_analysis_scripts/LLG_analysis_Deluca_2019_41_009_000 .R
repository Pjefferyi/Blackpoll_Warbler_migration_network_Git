require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(devtools)
library(remotes)
library(anytime)
library(ggmap)
library(FLightR)
library(LLmig)

#Sys.getenv("GITHUB_PAT")
#Sys.unsetenv("GITHUB_PAT")

#instal SGAT (dependency for TWGeos)
#devtools::install_github("SWotherspoon/SGAT")
library(SGAT)

#install the R package TwGeos from github to define twilight events 
#devtools::install_github("SLisovski/TwGeos")
library(TwGeos)

# TwGeos manual: https://rdrr.io/github/slisovski/TwGeos/

#other important packages 
#install_github("SLisovski/GeoLocTools")
library(GeoLocTools)
setupGeolocation()

# Extract data for a single bird in the data set from Deluca et al. 2019 =======
