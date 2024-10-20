# Script to assemble the  tracking data of all birds used in the geolocator analysis 
library(tidyverse)
library(TwGeos)
library(anytime)

# load reference data 
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/Geolocator_reference_data_consolidated.csv")

# Create a list of path to all files with location data
folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)
geo_names <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data")
names(folder_paths) <- geo_names

# filter to include only geolocators used in the Deluca et al. 2018 study, including study site that were ultimately note used in our own analysis
geos_prep <- ref.data %>% dplyr::filter(year(deploy.on.date) == 2016) %>% add_row(geo.id = c("4068_065", "4210_010", "4068_017"),
                                                                                  tag.model = c("MK-6","MK-6","MK-6"),
                                                                                  Standard.geo.id = c("4068-065", "4210-010", "4068-017"))
folder_paths <- folder_paths[which(names(folder_paths) %in% c(geos_prep$geo.id))]
geo_names <- geo_names[geo_names %in% c(geos_prep$geo.id)]

# create a table with all the light level measurements for each bird
twl_df <- data.frame()

for (i in seq(1:length(folder_paths))){

  # read twl data
  if (geos_prep[geos_prep$geo.id == geo_names[i],]$tag.model %in% c("MK-6",  "MK-6040",   "ML-6740")){
  lig <- readLig(paste0(folder_paths[i],"/Raw_light_data_", geo_names[i], ".lig"), skip = 1)

  lig <- lig %>% dplyr::select(Date, Light)

  } else {
  lig <- readMTlux(paste0(folder_paths[i],"/Raw_light_data_", geo_names[i], ".lux"), skip = 20)
  }

  # Add identifying info
  lig <- lig %>% mutate(geo_id = geos_prep[geos_prep$geo.id == geo_names[i],]$Standard.geo.id,  taxon = "setophaga striata", study_name = "A boreal songbird's migration across North America and the Atlantic Ocean, Setophaga striata",
                        sensor_type = "solar-geolocator")

  # Append to global dataset
  twl_df <- rbind(twl_df, lig)
}

# Save data for export
write_csv(twl_df, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/MoveBank_data_exports/Deluca_et_al_2019_geolocator_raw_data.csv")

# Create a list of path to all files with location data
folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)
geo_names <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data")
names(folder_paths) <- geo_names

# filter to include only geolocators deployed and recovered after 2017 
geos_prep <- ref.data %>% dplyr::filter(year(deploy.on.date) > 2016)

folder_paths <- folder_paths[which(names(folder_paths) %in% c(geos_prep$geo.id))]
geo_names <- geo_names[geo_names %in% c(geos_prep$geo.id)]

# create a table with all the light level measurements for each bird
twl_df <- data.frame()

for (i in seq(1:length(folder_paths))){

  # read twl data
  if (geos_prep[geos_prep$geo.id == geo_names[i],]$tag.model %in% c("MK-6",  "MK-6040",   "ML-6740")){
  lig <- readLig(paste0(folder_paths[i],"/Raw_light_data_", geo_names[i], ".lig"), skip = 1)

  lig <- lig %>% dplyr::select(Date, Light)

  } else {
  lig <- readMTlux(paste0(folder_paths[i],"/Raw_light_data_", geo_names[i], ".lux"), skip = 20)
  }

  # Add identifying info
  lig <- lig %>% mutate(geo_id = geo_names[i], "individual-taxon-canonical-name" = "setophaga striata", study_name = "Range-wide post- and pre-breeding migratory networks of a declining Neotropical-Nearctic migratory bird, the blackpoll warbler - Data collected after 2017",
                        "sensor-type" = "solar-geolocator") %>% dplyr::rename(timestamp = Date, "gls:light-level" = Light)
    
  # Append to global dataset
  twl_df <- rbind(twl_df, lig)
}

# Save data for export
write_csv(twl_df, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/MoveBank_data_exports/Post_2018_geolocator_raw_data.csv")
