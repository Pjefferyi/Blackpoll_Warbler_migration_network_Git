 findLocData <- function(geo.ids = c(), check_col_length = F, ref_path = NA, with_edits = c()){
  
  # Create a list of path to all files with location data 
  folder_paths <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", full.names = T)
  geo_names <- list.files("/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data")
  
  location_set <- data.frame()
  
  # Check the number of columns in the dataset for each geolocator 
  if (check_col_length == T){
    for (i in seq(1:length(folder_paths))){
      if (geo_names[i] %in% with_edits & geo_names[i] %in% geo.ids){
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary_fall_edit.csv")
        load(file = file_path)
        print(ncol(sm.fall.edit))
      }
      if (geo_names[i] %in% geo.ids){
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary.csv")
        load(file = file_path)
        print(folder_paths[i])
        print(ncol(sm))
        }
    }
    return()
  }
  
  
  for (i in seq(1:length(folder_paths))){
  # load the data from each file and add it to dataset if it is in geo_ids 
    if (geo_names[i] %in% geo.ids){
      # some of the data has edits during durign the fall transoceanic flight
      if (geo_names[i] %in% with_edits){
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary_fall_edit.csv")
        load(file = file_path)
        location_set <- rbind(location_set, sm.fall.edit)
      } else {
        file_path <- paste0(folder_paths[i], "/",geo_names[i],"_SGAT_GroupedThreshold_summary.csv")
        load(file = file_path)
        location_set <- rbind(location_set, sm)
      } 
    }
  }
  
  # Add reference data 
  if (!is.na(ref_path)){
    
    ref.data <- read.csv(ref_path)
    
    # Only retain relevant rows
    # ref.data <- ref.data[,c("geo.id",
    #                         "deploy.latitude",
    #                         "deploy.longitude",
    #                         "study.site",
    #                         "Range_region")]
    
    #Join the region and location data (inner join)
    location_set <- merge(location_set, ref.data, by.x = "geo_id", by.y = "geo.id")
  }
  
  return(location_set)
 }
 
 
 
 path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"   
 r3 <- findLocData(geo.ids = c("V8757_010", "V8296_004"), check_col_length = F, ref_path = path)