library(tidyverse)

# Load location database with all geolocators
geo.all <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/All.locations.csv")

# List of geolocator IDs used in the analysis
# Some geolocators in the data files were not used 
geo_names <- unique(geo.all$geo_id)

# Script to calculate summary statisics of the twilight editing process for all geolocators in the dataset 

# Create a list of path to all files with location data
folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)
folder_names <- list.dirs("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = F)

#create dataframe with folder paths and remove those for geolocators that were not used in their analysis 
folder_df <- data.frame(dir = folder_names[2:length(folder_names)], path = folder_paths) %>%
  dplyr::filter(dir %in% geo_names)

# dataframe that will contains twl data for all birds 
ag_twl <- data.frame()

for (i in seq(1:length(folder_df$dir))){
    file_path <- paste0(folder_df$path[i], "/",folder_df$dir[i],"_twl_times_edited.csv")
    twl_edited <- read.csv(file = file_path) %>% 
      mutate(geo_id = folder_df$dir[i])
    ag_twl <- ag_twl %>% rbind(twl_edited) %>%
      dplyr::select(geo_id, everything())
    print(folder_df$path[i])
}

# calculate the statistics of edits made using the twilight edit function
ag_delete_stat <- ag_twl %>% count(geo_id, Deleted) %>% complete(geo_id, Deleted, fill = list(n = 0)) %>%
  arrange(geo_id) %>% pivot_wider(names_from = Deleted, values_from = n) %>% mutate(prop_deleted = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  summarize(mean_deleted = mean(prop_deleted), sd_deleted = sd(prop_deleted), n_deleted = n(),
            min_deleted = min(prop_deleted), max_deleted = max(prop_deleted))

ag_edit_stat <- ag_twl %>% count(geo_id, Edited) %>% complete(geo_id, Edited, fill = list(n = 0)) %>%
  arrange(geo_id) %>% pivot_wider(names_from = Edited, values_from = n) %>% mutate(prop_edited = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  summarize(mean_editd = mean(prop_edited), sd_editd = sd(prop_edited), n_edited = n(),
            min_editd = min(prop_edited), max_editd = max(prop_edited))
