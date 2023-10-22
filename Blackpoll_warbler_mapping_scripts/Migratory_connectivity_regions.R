# Plot of the breeding and nonnonbreeding regions used to calculated the MC metric

# Load geolocator analysis helper functions 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# run the network preparation script to get the nonbreeding locations for individuals in the fall network
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Network_construction.R")

# Plot of breeding regions ####

# fall network breeding points
fall.br <- fall.breed %>% arrange(geo_id)

# Fall network breeding regions 
fall.br.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_abundance_regions_adjusted.shp") %>%
  st_transform(crs(wrld_simpl)) %>%
  st_cast("MULTIPOLYGON")

# Create the plot in base R
plot(as_Spatial(fall.br.regions), lwd = 0.5, xlim = c(-170, -25), col = adjustcolor(c("darkgreen", "darkblue", "purple"), alpha = 0.5), ylim = c(50, 70))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], lwd = 1, add = T)

points(fall.br$Lon.50., fall.br$Lat.50., pch = 16, cex = 2, col = "white")
points(fall.br$Lon.50., fall.br$Lat.50., pch = 1, cex = 2, col = "black")


# Create the plot in ggplot
ggplot(st_as_sf(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),]))+
  geom_sf(colour = "black", fill = "lightgray") +
  geom_sf(data = fall.br.regions, aes(fill = region), alpha = 0.8)+
  scale_fill_discrete(name = "Breeding \n regions", label = c("Eastern", "Central", "Western"))+
  coord_sf(xlim = c(-170, -35), ylim = c(40, 75)) +
  geom_point(data = fall.br, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), col = "white", fill = "black", pch = 21) +
  theme_bw() +
  theme(text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab("Longitude")+
  ylab("Latitude")


# Plot of nonbreeding regions ####

# Load the nonbreeding regions  
fall.nbr.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Migratory connectivity_regions/Data/NonbreedingregionsV3.shp") %>%
  st_transform(crs(wrld_simpl))%>%
  st_cast("MULTIPOLYGON")

# Get the first nonbreeding site of every bird 
fall.nbr <- fall.stat.ab %>% group_by(geo_id) %>% 
  filter(sitenum == max(sitenum), !is.na(StartTime), !is.na(EndTime)) %>%
  arrange(geo_id) 

# Create the plot 
plot(as_Spatial(fall.nbr.regions), col = adjustcolor(c("lightblue", "orange"), alpha = 0.8), lwd = 0.5,
     xlim = c(-90, -30), ylim = c(-10, 15))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     col = NA, lwd = 1, add = T)

points(fall.nbr$Lon.50., fall.nbr$Lat.50., pch = 16, cex = 2, col = "white")
points(fall.nbr$Lon.50., fall.nbr$Lat.50., pch = 1, cex = 2, col = "black")


# Create the plot in ggplot
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "lightgray") +
  geom_sf(data = fall.nbr.regions, aes(fill = region), alpha = 0.8)+
  scale_fill_discrete(name = "Nonbreeding \n regions")+
  coord_sf( xlim = c(-90, -30), ylim = c(-10, 15)) +
  geom_point(data = fall.nbr, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), col = "white", fill = "black", pch = 21) +
  theme_bw() +
  theme(text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab("Longitude")+
  ylab("Latitude")
