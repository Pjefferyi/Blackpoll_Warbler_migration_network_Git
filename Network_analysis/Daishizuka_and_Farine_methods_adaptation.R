


# Load fall stat and cluster data 
fall.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.data.csv")

smps <- 1000 # number of boostrap resamples

list 

for (i %in% seq(1:smps)){
  
  individuals <- sample(unique(fall.stat$geo_id), 
                        size = length(unique(fall.stat$geo_id)),
                        replace = T)
  
  new.data <- fall.stat %>% filter(geo_id %in% individuals)
  
  #create the new network
  # Add a column with inter-cluster connections to our location dataset
  new.data <- new.data%>% mutate(next.cluster = case_when(
    lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
    .default = NA))
  
  # calculate abundance propagation units for each individual 
  
  
  # Create new edge data 
  new.data.edge <- new.data %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                              Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                              Lat.2.5., Lat.97.5., EndTime,
                                              sitenum, duration, period,study.site,
                                              Range_region, NB_count, period,
                                              site_type) %>%
  filter(!is.na(next.cluster))
  
  new.vertices <- data.frame("vertex" = seq(1, max(new.data.edge$cluster)))
  
  new.graph <- graph_from_data_frame(new.data.edge, directed = T, vertices = new.vertices)
  
  

}
