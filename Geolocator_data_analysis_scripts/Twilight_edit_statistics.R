library(tidyverse)

# Script to calculate summary statisics of the twilight editing process for all geolocators in the dataset 

# Create a list of path to all files with location data
folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)

# List of geo names
geo_names <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data")

# dataframe that will contains twl data for all birds 
ag_twl <- data.frame()

for (i in seq(1:length(folder_paths))){
    file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary.csv")
    load(file = file_path)
    print(folder_paths[i])
    print(ncol(sm))
}