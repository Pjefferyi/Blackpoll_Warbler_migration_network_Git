################################################################################
# Relation between the breeding and nonbreeding location of blackpoll warblers
################################################################################

# load libraries 
library(ggplot2)
library(anytime)
library(performance)
library(ggpubr)

# First extract analyzied geolocator data 
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
                                   "blpw09",
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
                                   #"E",
                                   "D"), check_col_length = F)

# Get reference data 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# Obtain breeding location for each geolocator 
geo.br.loc <- geo.all %>% group_by(geo_id) %>%
  summarize(breed.lat = unique(deploy.latitude ), breed.lon = unique(deploy.longitude))

# for WRMA04173 (tracked from nonbreeding ground) we will use the arrival location in NA
wrma.x <- geo.all[(geo.all$geo_id == "WRMA04173"),] %>% filter(sitenum == max(sitenum)) %>% dplyr::select(Lon.50.)
wrma.y <- geo.all[(geo.all$geo_id == "WRMA04173"),] %>% filter(sitenum == max(sitenum)) %>% dplyr::select(Lat.50.)

geo.br.loc[(geo.br.loc$geo_id == "WRMA04173"),]$breed.lat <- wrma.y$Lat.50. 
geo.br.loc[(geo.br.loc$geo_id == "WRMA04173"),]$breed.lon <- wrma.x$Lon.50.

# Obtain first nonbreeding location for each geolocator
geo.nbr1.loc <- geo.all %>% dplyr::filter(period == "Non-breeding period", sitenum > 0, duration > 14) %>% group_by(geo_id) %>%
  dplyr::filter(sitenum == min(sitenum)) %>%
  dplyr::select(geo_id, nbr1.lat = Lat.50., nbr1.lon = Lon.50., nbr.duration = duration)

# Obtain last nonbreeding location for each geolocator
geo.nbr2.loc <- geo.all %>% dplyr::filter(period == "Non-breeding period", sitenum > 0, duration > 14) %>% group_by(geo_id) %>%
  dplyr::filter(sitenum == max(sitenum)) %>%
  dplyr::select(geo_id, nbr2.lat = Lat.50., nbr2.lon = Lon.50., nbr.duration = duration)

# merge the datasets
loc_list <- list(geo.br.loc, geo.nbr1.loc , geo.nbr2.loc)
geo.loc.data <- loc_list %>% reduce(full_join, by='geo_id') %>% 
  merge(ref_data, by.x = "geo_id", by.y = "geo.id")

# Plot the first nonbreeding longitude against the breeding longitude 
long.plot1 <- ggplot(data = geo.loc.data, aes(x = breed.lon, y = nbr1.lon)) + 
  geom_point(data = geo.loc.data, aes(x = breed.lon, y = nbr1.lon)) + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "First nonbreeding site longitude") # +
  #geom_text(label = geo.loc.data$geo_id, nudge_x = 5, nudge_y = 1)

# run a regression model
mod1 <- lm(nbr1.lon ~ breed.lon, data = geo.loc.data)
summary(mod1)
#check_model(mod1)

# Plot the second nonbreeding longitude against the breeding longitude 
long.plot2 <- ggplot(data = geo.loc.data, aes(x = breed.lon, y = nbr2.lon)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme(text = element_text(size=14)) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "second nonbreeding site longitude")

# run a regression model
mod2 <- lm(nbr2.lon ~ breed.lon, data = geo.loc.data)
summary(mod2)
#check_model(mod2)

# Plot the first nonbreeding latitude against the breeding latitude
lat.plot1 <- ggplot(data = geo.loc.data, aes(x = breed.lat, y = nbr1.lat)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "First nonbreeding site latitude")

# run a regression model
mod3 <- lm(nbr1.lat ~ breed.lat, data = geo.loc.data)
summary(mod3)
#check_model(mod3)

# Plot the first nonbreeding latitude against the breeding latitude
lat.plot2 <- ggplot(data = geo.loc.data, aes(x = breed.lat, y = nbr2.lat)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "Second nonbreeding site latitude")

# run a regression model
mod4 <- lm(nbr2.lat ~ breed.lat, data = geo.loc.data)
summary(mod4)
#check_model(mod4)

# Combine all plots into one 
ggarrange(long.plot1, long.plot2)
ggarrange(lat.plot1, lat.plot2)
