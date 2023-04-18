#Check that the duration of the tracking period for all geolocators

# retrieve path to all raw light data files 

# list of paths 
list <-list.files("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", pattern = "\\.lig", recursive = T, 
                  full.names = T)

#list of filenames 
list.short <- list.files("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geolocator_data", pattern = "\\.lig", recursive = T) 

# Create a dataframe to store the start and end dates fo the tracking periods 

dur.data <- data.frame(File_name = character(), 
                       Start_date = as.Date(x = integer(0), origin = "1970-01-01"),
                       End_date = as.Date(x = integer(0), origin = "1970-01-01"))

for (e in list){
  
  lig <- readLig(e, skip = 1)
  
  dur.data[nrow(dur.data) + 1,] <- list(list.short[nrow(dur.data) + 1], 
                                     lig$Date[1],
                                     lig$Date[nrow(lig)])
}

# display the results 
dur.data
