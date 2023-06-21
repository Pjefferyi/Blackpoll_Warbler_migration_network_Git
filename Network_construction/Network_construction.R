#Load libraries
library(tidyverse)
library(igraph)
library(lubridate)
library(anytime)
library(ggplot2)

# Clustering breeding and nonbreeding location #################################

# Either import a file with manual clusters, or create cluster in R


# Create clusters in R 
geo.all <- findLocData(geo.ids = c("V8757_010",
                                   "V8296_004",
                                   "V8296_005",
                                   "V8296_006",
                                   "V8757_055",
                                   "V8757_018",
                                   "V8757_021",
                                   "V8296_015",
                                   "V8296_017",
                                   "V8296_026",
                                   "V8296_025",
                                   "V8296_007",
                                   "V8296_008",
                                   "V8757_019",
                                   "V8757_096",
                                   "V8757_134",
                                   "V8757_029",
                                   "V8757_078",
                                   "blw09",
                                   "blpw12",
                                   "3254_001",
                                   "4068_014",
                                   "blpw14",
                                   "3254_003",
                                   "3254_008",
                                   "3254_011",
                                   "3254_057",
                                   "blpw15",
                                   "blpw25",
                                   "4105_008",
                                   "4105_009",
                                   "4105_016",
                                   "4105_017",
                                   "4210_002",
                                   "4210_004",
                                   "4210_006",
                                   "4210_010",
                                   "WRMA04173",
                                   "A",
                                   "B",
                                   "C",
                                   "D",
                                   "E"), check_col_length = F, ref_path = ref_path)


# Generate edge list for the network ###########################################

Cl.data 