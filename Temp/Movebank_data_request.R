# This script is to build the dataset requested for addition to the Audubon Migration explorer 
library(tidyverse)
library(lubridate)
library(anytime)
library(rnaturalearth)
library(terra)
library(sf)

# Create a list of path to all files with location data 
folder_paths <- list.files("C:/Users/Jelan/Github/Geolocator_data", full.names = T)
geo_names <- list.files("C:/Users/Jelan/Github/Geolocator_data")

#load the geolocator reference data/metadata
ref.data <- read.csv("C:/Users/Jelan/Github/Blackpoll_Warbler_migration_network_Git/Analysis_output_data/Geolocator_reference_data_consolidated.csv")

#Vector of geolocators for which data will be extracted
geo.ids <- c("V8757_010",
             "V8296_004",
             "V8296_005",
             "V8296_006",
             "V8757_055",
             "V8757_018",
             "V8296_015",
             "V8296_017",
             "V8296_021",
             "V8296_026",
             "V8296_025",
             "V8296_007",
             "V8296_008",
             "V8757_019",
             "V8757_096",
             "V8757_134",
             "V8757_029",
             "V8757_078",
             "V7638_001",
             "V7638_005",
             "V7638_009",
             "V7638_010",
             "V7638_011",
             "WRMA04173")


#Empty dataframe that will be used to store location data
location_set <- data.frame()

# iterate through the folders to obtain the necessary data 
for (i in seq(1:length(folder_paths))){
  # load the data from each file and add it to dataset if it is in geo_ids 
  if (geo_names[i] %in% geo.ids){
    
    if (ref.data[ref.data$geo.id == geo_names[i],]$Clock_drift_edits == F){
      print(geo_names[i])
      #extract threshold location data 
      load(paste0(folder_paths[i], "/",geo_names[i],"_initial_path.csv"), verbose = T)
      x0 <- as.data.frame(x0)
      
      # Add column with the geo ID
      x0$geo_id <- geo_names[i]
      
      # Add times 
      twl <- read.csv(paste0(folder_paths[i],"/", geo_names[i], "_twl_times_edited.csv"))
      x0$twilight <- twl$Twilight
      x0$Rise <- twl$Rise
    
      } else {
        print(geo_names[i])
        #extract threshold location data 
        load(paste0(folder_paths[i], "/",geo_names[i],"_initial_path__clock_drift_adjustment.csv"), verbose = T)
        x0 <- as.data.frame(x0)
        # Add column with the geo ID
        x0$geo_id <- geo_names[i]
        
        # Add times 
        twl <- read.csv(paste0(folder_paths[i],"/", geo_names[i], "_twl_times_edited.csv"))
        x0$twilight <- twl$Twilight
        x0$Rise <- twl$Rise
      }
    location_set <- rbind(location_set, x0)
    }
}

# Add reference data 
location_set <- merge(location_set, ref.data[, c("geo.id", "animal.taxon", "band.id", "Age.status", "animal.sex", "study.site")], by.x = "geo_id", by.y = "geo.id")

# Some geolocators have errors  in the year of deployment
unique(filter(location_set[year(location_set$twilight) %in% c(2012, 2013),])$geo_id)

# Fix these differences
location_set$twilight <- anytime(location_set$twilight)
year(location_set[year(location_set$twilight) %in% c(2012, 2013),]$twilight) <- year(location_set[year(location_set$twilight) %in% c(2012, 2013),]$twilight) + 7

# Save
write.csv(location_set, "C:/Users/Jelan/OneDrive/Desktop/Blpw_location_estimates.csv")
save(location_set, file = "C:/Users/Jelan/OneDrive/Desktop/Blpw_location_estimates.R")

# Plot a sample 
world_countries <- ne_countries(type = "countries", scale = "large")
America <- world_countries[((world_countries$continent %in%  c("North America", "South America") | world_countries$admin == "France") & world_countries$admin != "Greenland"),]

eg <- location_set[location_set$geo_id ==  "V8296_005", ]
plot(as_Spatial(America))
points(eg$lon, eg$lat)
