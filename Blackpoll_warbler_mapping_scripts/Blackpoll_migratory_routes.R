# Load helper functions & packages 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis/Geolocator_analysis_helper_functions.R")

# Plot all geolocator tracks for the spring ####################################

geo.list <- c("V8757_010",
              "V8296_004",
              "V8296_005",
              "V8296_006",
              "V8757_055",
              "V8757_018",
              #"V8757_021",
              #"V8296_015",
              "V8296_017",
              "V8296_026",
              #"V8296_025",
              #"V8296_007",
              #"V8296_008",
              "V8757_019",
              "V8757_096",
              "V8757_134",
              "V8757_029",
              "V8757_078",
              "blpw12",
              "3254_001",
              "4068_014",
              "blpw14",
              "3254_003",
              "3254_008",
              "blpw15")
  

# retrive location data 
geo.all <- findLocData(geo.ids = geo.list , check_col_length = F)

# retrieve region data (Western and central, or Eastern)
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

region.data <- ref.data[,c("geo.id", "Range_region")]

#Join the region and location data (inner join)
geo.all <- merge(geo.all, region.data, by.x = "geo_id", by.y = "geo.id")

# extract breeding location 
breed.site.data <- ref.data[(ref.data$geo.id %in% geo.list),]

# extract spring routes 
geo.spring.all <- geo.all[(geo.all$period %in% c("Pre-breeding migration", "Non-breeding period")),]

# Now we can plot the spring migration routes 
ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
  geom_sf() +
  coord_sf() +
  geom_path(data = geo.spring.all, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, color = Range_region), linewidth = 0.3) +
  #geom_errorbar(data = data, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") + 
  #geom_errorbar(data = data, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") + 
  geom_point(data = breed.site.data, mapping = aes(x = deploy.longitude, y = deploy.latitude), size = 1.4, col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 12))

# plot of spring stationary locations
geo.spring.stat <- geo.spring.all[(geo.spring.all$sitenum > 0),]

ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = geo.spring.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "black", alpha = 0.5) + 
  geom_errorbar(data = geo.spring.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "black", alpha = 0.5) + 
  geom_path(data = geo.spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, color = Range_region), linewidth = 0.3, alpha = 0.8) +
  geom_point(data = geo.spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, color = Range_region, fill = Range_region), size = 1.3, shape = 23) +
  geom_point(data = breed.site.data, mapping = aes(x = deploy.longitude, y = deploy.latitude), size = 1.4, col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 14))

# plot nonbreeding areas ######################################################

geo.list <- c("V8757_010",
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
              "blpw12",
              "3254_001",
              "4068_014",
              "blpw14",
              "3254_003",
              "3254_008",
              "blpw15")

# retrieve location data 
geo.all <- findLocData(geo.ids = geo.list , check_col_length = F)

# extract nonbreeding locations
geo.nbr.all <- geo.all[(geo.all$period %in% c("Non-breeding period")),]
geo.nbr.stat<- geo.nbr.all[(geo.nbr.all$sitenum > 0),]

# crop map of south america
crop.nbr <- st_crop(spData::world[(spData::world$continent %in% c( "South America")),],
                    xmin = -90, xmax = -30, ymin = -20, ymax = 20)

# Now we can plot the nonbreeding locations
ggplot(crop.nbr) +
  geom_sf() +
  coord_sf() +
  geom_path(data = geo.nbr.stat, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id), linewidth = 0.3) +
  geom_point(data = geo.nbr.stat, mapping = aes(x = Lon.50., y = Lat.50., colour= geo_id, fill = geo_id), size = 1.3, shape = 23) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 12))


# Plot all geolocator tracks for the fall  #####################################

geo.list <- c("V8757_010",
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
              "blpw12",
              "3254_001",
              "4068_014",
              "blpw14",
              "3254_003",
              "3254_008",
              "blpw15")

# retrive location data 
geo.all <- findLocData(geo.ids = geo.list , check_col_length = F)

# retrieve region data (Western and central, or Eastern)
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

region.data <- ref.data[,c("geo.id", "Range_region")]

#Join the region and location data (inner join)
geo.all <- merge(geo.all, region.data, by.x = "geo_id", by.y = "geo.id")

# extract breeding location 
breed.site.data <- ref.data[(ref.data$geo.id %in% geo.list),]

# extract fall routes 
geo.fall.all <- geo.all[(geo.all$period %in% c("Post-breeding migration", "Non-breeding period")),]

# Now we can plot the fall migration routes 
ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
  geom_sf() +
  coord_sf() +
  geom_path(data = geo.fall.all, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, color = Range_region), linewidth = 0.3) +
  #geom_errorbar(data = data, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") + 
  #geom_errorbar(data = data, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") + 
  geom_point(data = breed.site.data, mapping = aes(x = deploy.longitude, y = deploy.latitude), size = 1.4, col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 12))

# plot of fall stationary locations
geo.fall.stat <- geo.fall.all[(geo.fall.all$sitenum > 0),]

ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),]) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = geo.fall.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "black", alpha = 0.5) + 
  geom_errorbar(data = geo.fall.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "black", alpha = 0.5) + 
  geom_path(data = geo.fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, color = Range_region), linewidth = 0.3, alpha = 0.8) +
  geom_point(data = geo.fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, color = Range_region, fill = Range_region), size = 1.3, shape = 23) +
  geom_point(data = breed.site.data, mapping = aes(x = deploy.longitude, y = deploy.latitude), size = 1.4, col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 14))

# Plot migratory routes passing over the carribeans ############################

geo.list <- c("V8296_004",
              "V8296_005",
              "V8296_006",
              "V8757_018",
              "V8757_021",
              "V8296_015",
              "V8296_017",
              "V8296_026",
              "V8296_025",
              "V8296_007",
              "V8296_008",
              "V8757_096",
              "V8757_078",
              "V8757_134")

# retrive location data 
geo.all <- findLocData(geo.ids = geo.list , check_col_length = F)

# extract fall routes and stat locations
geo.fall.all <- geo.all[(geo.all$period %in% c("Post-breeding migration", "Non-breeding period")),]
geo.fall.stat <- geo.fall.all[(geo.fall.all$sitenum > 0),]

# crop map to show carribean crossing
crop.carib <- st_crop(spData::world[(spData::world$continent %in% c("North America", "South America")),],
                    xmin = -90, xmax = -45, ymin = -10, ymax = 50)

# retrieve region data (Western and central, or Eastern)
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

# extract breeding location 
breed.site.data <- ref.data[(ref.data$geo.id %in% geo.list),]

ggplot(crop.carib) +
  geom_sf() +
  coord_sf() +
  geom_path(data = geo.fall.stat, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id), linewidth = 0.3, alpha = 0.8) +
  geom_point(data = geo.fall.stat, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id, fill = geo_id), size = 1.3, shape = 23) +
  geom_point(data = breed.site.data, mapping = aes(x = deploy.longitude, y = deploy.latitude), size = 1.4, col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #theme(legend.position = "none") +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 14))

# Plot migratory routes passing over the Southeast America ############################

geo.list <- c("V8757_019", "V8757_029", "4068_014", "blpw14")

# retrive location data 
geo.all <- findLocData(geo.ids = geo.list , check_col_length = F)

# extract fall routes and stat locations
geo.S.US.all <- geo.all[(geo.all$period %in% c("Post-breeding migration", "Non-breeding period")),]
geo.S.US.stat <- geo.S.US.all[(geo.S.US.all$sitenum > 0),]

# crop map to show carribean crossing
crop.S.US <- st_crop(spData::world[(spData::world$continent %in% c("North America", "South America")),],
                      xmin = -170, xmax = -50, ymin = -10, ymax = 70)

# retrieve region data (Western and central, or Eastern)
ref.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv")

# extract breeding location 
breed.site.data <- ref.data[(ref.data$geo.id %in% geo.list),]

ggplot(crop.S.US) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = geo.S.US.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "black", alpha = 0.5) + 
  geom_errorbar(data = geo.S.US.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "black", alpha = 0.5) + 
  geom_path(data = geo.S.US.stat, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id), linewidth = 0.3, alpha = 0.8) +
  geom_point(data = geo.S.US.stat, mapping = aes(x = Lon.50., y = Lat.50., color = geo_id, fill = geo_id), size = 1.3, shape = 23) +
  geom_point(data = breed.site.data, mapping = aes(x = deploy.longitude, y = deploy.latitude), size = 1.4, col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #theme(legend.position = "none") +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 14))

# Empty map of America #########################################################

ggplot(spData::world[(spData::world$continent %in% c("North America", "South America")),], fill = ) +
  geom_sf(colour = "black", fill = "white", linewidth=0.3) +
  coord_sf() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
