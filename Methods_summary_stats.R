# Thesis methods details 

# load necessary libraries 
library(performance)
library(DHARMa)
library(MASS)

# Load the helper functions script  
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Load geolocator reference data 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# location data (Only for geolocator used in the analysis )
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
              "V7638_001",
              "V7638_005",
              "V7638_009",
              "V7638_010",
              "V7638_011",
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
              "D")

geo.all <- findLocData(geo.ids = geo.list, check_col_length = F)

# filter reference data to include only geolocators used in the analysis 
analysis_ref <- ref_data %>% dplyr::filter(geo.id %in% unique(geo.all$geo_id))

# mean Zenith angle measured using in-habitat calibration ----
IH_zentih <- mean(analysis_ref$In_habitat_median_zenith_angle, na.rm = T)

# mean zenith angle measured using Hill-Ekstrom Calibration ----
He_zenith <- mean(analysis_ref$Hill_Ekstrom_median_angle, na.rm = T)

# mean difference between the Ih and HE zenith ----
zenith_diff <- mean(analysis_ref$In_habitat_median_zenith_angle - analysis_ref$Hill_Ekstrom_median_angle, na.rm = T)
hist(analysis_ref$In_habitat_median_zenith_angle - analysis_ref$Hill_Ekstrom_median_angle, breaks = 20)

# Linear interpolation period for latitude in preliminary location estimates 

# Get the solar declination angle for all days in a year 
doy <- seq(anytime("2023-06-01"), anytime("2023-12-31"), by = "day")
solar.angle  <- solar(doy)$sinSolarDec

# Calculate  the number of days censored prior to and after the equinox based on the value of the tol parameter 
tol_values <- seq(0, 0.25, by = 0.001) 
fall_equi <- anytime("2023-09-22")
cens_data <- data.frame(tol_values)

cens_data <- cens_data %>% 
  rowwise() %>%
  mutate(days_cens = length(solar.angle[solar.angle > 0 - tol_values & solar.angle < 0 + tol_values])/2)

# calculate the number of days censored for each geolocator 
analysis_ref <- analysis_ref %>% rowwise() %>% mutate(tol_days = list(cens_data[tol_values == tol,]$days_cens),
                                                      tol_days = tol_days[[1]])
# Range of days 
range(analysis_ref$tol_days)

# Mean and SD
mean(analysis_ref$tol_days)
sd(analysis_ref$tol_days)

# Results #####################################################################

# Load fall and spring stationary sites 
fall.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.data.csv")
spring.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.stationary.data.csv")

## Number of stationary locations ----

# Number of fall stopovers used before arrival at first nonbreeding site
fall.num <- fall.stat %>% group_by(geo_id) %>% summarise(num_stopovers = sum(site_type == "Stopover"))

# Number spring site used after departure from last nonbreeding site 
spring.num <- spring.stat %>% group_by(geo_id) %>% summarise(num_stopovers = sum(site_type == "Stopover"))

## Number migration distance ----

# Fall migration distance 
fall.dist <- fall.stat %>% group_by(geo_id) %>% mutate(Lon.50.next = lead(Lon.50.),
                                                      Lat.50.next = lead(Lat.50.)) %>%
  rowwise() %>%
  mutate(dists = distHaversine(c(Lon.50.,Lat.50.), c(Lon.50.next, Lat.50.next))) %>%
  ungroup() %>%
  group_by(geo_id) %>%
  filter(!is.na(dists)) %>%
  summarize(distance = sum(dists)/1000)

# Spring migration distance
spring.dist <- spring.stat %>% group_by(geo_id) %>% mutate(Lon.50.next = lead(Lon.50.),
                                                       Lat.50.next = lead(Lat.50.)) %>%
  rowwise() %>%
  mutate(dists = distHaversine(c(Lon.50.,Lat.50.), c(Lon.50.next, Lat.50.next))) %>%
  ungroup() %>%
  group_by(geo_id) %>%
  filter(!is.na(dists)) %>%
  summarize(distance = sum(dists)/1000)

## Number of nodes used ----

# Load fall and spring intra node movement datasets
fall.edge.df.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.intra.cluster.movements.csv")
spring.edge.df.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.intra.cluster.movements.csv")

# add node type info 
meta.fall.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")
meta.spring.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")

fall.edge.df.ab <- merge(fall.edge.df.ab, meta.fall.ab[,c("vertex","node.type")], by.x = "cluster", by.y = "vertex")
spring.edge.df.ab <- merge(spring.edge.df.ab, meta.spring.ab[,c("vertex","node.type")], by.x = "cluster", by.y = "vertex")

# number of nodes used in the fall network
fall.node.used.df <- fall.edge.df.ab %>% group_by(geo_id) %>% summarize(fall.nodes.occupied = length(unique(cluster, next.cluster)),
                                                                        fall.stopover.nodes.occupied = length(unique(cluster[node.type == "Stopover"], next.cluster[node.type == "Stopover"])),
                                                                        fall.nbr.nodes.occupied = length(unique(cluster[node.type == "Nonbreeding"], next.cluster[node.type == "Nonbreeding"])))

# number of nodes used in the spring network 
spring.node.used.df <- spring.edge.df.ab %>% group_by(geo_id) %>% summarize(spring.nodes.occupied = length(unique(cluster, next.cluster)),
                                                                        spring.stopover.nodes.occupied = length(unique(cluster[node.type == "Stopover"], next.cluster[node.type == "Stopover"])),
                                                                        spring.nbr.nodes.occupied = length(unique(cluster[node.type == "Nonbreeding"], next.cluster[node.type == "Nonbreeding"])))

# merge the information collected with the reference dataset 
analysis_ref <- merge(analysis_ref, fall.node.used.df, by.x = "geo.id", by.y = "geo_id", all = T)
analysis_ref <- merge(analysis_ref, spring.node.used.df, by.x = "geo.id", by.y = "geo_id", all = T)

## Statistics of node usage ----

# Average number of nodes used in fall network
range(analysis_ref$fall.nodes.occupied, na.rm =  T)
mean(analysis_ref$fall.nodes.occupied, na.rm =  T)
se <- sd(analysis_ref$fall.nodes.occupied, na.rm =  T)/sqrt(length(analysis_ref$fall.nodes.occupied[!is.na(analysis_ref$fall.nodes.occupied)]))

# Average number of nodes used in the spring network 
range(analysis_ref$spring.nodes.occupied, na.rm =  T)
mean(analysis_ref$spring.nodes.occupied, na.rm =  T)
se <- sd(analysis_ref$spring.nodes.occupied, na.rm =  T)/sqrt(length(analysis_ref$spring.nodes.occupied[!is.na(analysis_ref$spring.nodes.occupied)]))

# correlation between number of nodes used and longitude of breeding site 
mod.fall.node <- glm(fall.nodes.occupied ~ deploy.longitude, data = analysis_ref, family = gaussian(link = "identity"))
plot(fall.nodes.occupied ~ deploy.longitude, data = analysis_ref)
summary(mod.fall.node)
check_model(mod.fall.node)

mod.spring.node <- glm(spring.nodes.occupied ~ (deploy.longitude), data = analysis_ref, family = gaussian(link = "identity"))
plot(spring.nodes.occupied ~ deploy.longitude, data = analysis_ref)
summary(mod.spring.node)
check_model(mod.spring.node)

# test assessing the usage of more than one nonbreeding site (whether individuals made stopovers in the nonbreeding range)
analysis_ref <- analysis_ref %>% mutate(fall.nbr.stopover = ifelse(fall.nbr.nodes.occupied >1, 1, 0),
                                        spring.nbr.stopover = ifelse(spring.nbr.nodes.occupied >1, 1, 0))

nbr.stopover.mod.fall <- glm(fall.nbr.stopover ~ deploy.longitude, data = analysis_ref, binomial(link = "logit"))
boxplot(analysis_ref$deploy.longitude ~ as.factor(analysis_ref$fall.nbr.stopover))
summary(nbr.stopover.mod.fall)
check_model(nbr.stopover.mod.fall)

nbr.stopover.mod.spring <- glm(spring.nbr.stopover ~ deploy.longitude, data = analysis_ref, family = binomial(link = "logit"))
boxplot(analysis_ref$deploy.longitude ~ as.factor(analysis_ref$spring.nbr.stopover))
summary(nbr.stopover.mod.spring)
check_model(nbr.stopover.mod.spring)

simulationOutput <- simulateResiduals(fittedModel = nbr.stopover.mod.spring, plot = F)
plot(simulationOutput)

## Linear regression between the breeding and nonbreeding longitudes and latitudes -----

# for WRMA04173 (tracked from nonbreeding ground) we will use the arrival location in NA
wrma.x <- geo.all[(geo.all$geo_id == "WRMA04173"),] %>% filter(sitenum == max(sitenum)) %>% dplyr::select(Lon.50.)
wrma.y <- geo.all[(geo.all$geo_id == "WRMA04173"),] %>% filter(sitenum == max(sitenum)) %>% dplyr::select(Lat.50.)

analysis_ref[(analysis_ref$geo.id == "WRMA04173"),]$deploy.latitude <- wrma.y$Lat.50. 
analysis_ref[(analysis_ref$geo.id == "WRMA04173"),]$deploy.longitude <- wrma.x$Lon.50.

# Obtain first nonbreeding location for each geolocator
geo.nbr1.loc <- geo.all %>% dplyr::filter(period == "Non-breeding period", sitenum > 0, duration > 14, geo_id %in% unique(fall.stat$geo_id)) %>% group_by(geo_id) %>%
  dplyr::filter(sitenum == min(sitenum)) %>%
  dplyr::select(geo_id, nbr1.lat = Lat.50., nbr1.lon = Lon.50.)

# Obtain last nonbreeding location for each geolocator
geo.nbr2.loc <- geo.all %>% dplyr::filter(period == "Non-breeding period", sitenum > 0, duration > 14, geo_id %in% unique(spring.stat$geo_id)) %>% group_by(geo_id) %>%
  dplyr::filter(sitenum == max(sitenum)) %>%
  dplyr::select(geo_id, nbr2.lat = Lat.50., nbr2.lon = Lon.50.)

# merge with the dataset of reference data 
loc_list <- list(geo.nbr1.loc , geo.nbr2.loc)
analysis_ref <- loc_list %>% reduce(full_join, by='geo_id') %>% 
  merge(analysis_ref, by.x = "geo_id", by.y = "geo.id", all = T)

# Plot the first nonbreeding longitude against the breeding longitude 
long.plot1 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr1.lon)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "First nonbreeding site longitude") # +
#geom_text(label = geo.loc.data$geo_id, nudge_x = 5, nudge_y = 1)

# run a regression model
mod1 <- lm(nbr1.lon ~ deploy.longitude, data = analysis_ref)
summary(mod1)
check_model(mod1)

# Plot the second nonbreeding longitude against the breeding longitude 
long.plot2 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr2.lon)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme(text = element_text(size=14)) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "second nonbreeding site longitude")

# run a regression model
mod2 <- lm(nbr2.lon ~ deploy.longitude, data = analysis_ref)
summary(mod2)
check_model(mod2)

# Plot the first nonbreeding latitude against the breeding latitude
lat.plot1 <- ggplot(data = analysis_ref, aes(x = deploy.latitude, y = nbr1.lat)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "First nonbreeding site latitude")

# run a regression model
mod3 <- lm(nbr1.lat ~ deploy.latitude, data = analysis_ref)
plot(nbr1.lat ~ deploy.latitude, data = analysis_ref)
summary(mod3)
check_model(mod3)

# Plot the first nonbreeding latitude against the breeding latitude
lat.plot2 <- ggplot(data = analysis_ref, aes(x = deploy.latitude, y = nbr2.lat)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "Second nonbreeding site latitude")

# run a regression model
mod4 <- lm(nbr2.lat ~ deploy.latitude, data = analysis_ref)
plot(nbr2.lat ~ deploy.latitude, data = analysis_ref)
summary(mod4)
check_model(mod4)

x <- analysis_ref$nbr2.lon[!is.na(analysis_ref$nbr2.lon)]
y <- analysis_ref$deploy.longitude[!is.na(analysis_ref$nbr2.lon)]
cor(x,y, method = "pearson")

## network metric scores by nodes ----

#load network data
fall.gdata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/fall.graph.data.csv")
spring.gdata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/spring.graph.data.csv")

# merge with fall graph data with fall.stat region names 
fall.stat.regions <- fall.stat %>% mutate(cluster.region = ifelse(is.na(cluster.region), as.character(cluster), cluster.region)) %>%
  group_by(cluster) %>% summarise(cluster.region = unique(cluster.region))
fall.gdata <-  merge(fall.gdata, fall.stat.regions, by = "cluster")

# Plot betweenness scores
ggplot(data = fall.gdata, aes(y = betweenness, x = cluster.region, fill = node.type))+
  geom_col()+
  coord_flip()
  
# Plot bridge betweenness scores
ggplot(data = fall.gdata, aes(y = bridge.strength, x = cluster.region, fill = node.type))+
  geom_col()+
  coord_flip()

# There are no spring region names, but I may add some at a later date
# spring.stat.regions <- spring.stat %>% mutate(cluster.region = ifelse(is.na(cluster.region), as.character(cluster), cluster.region)) %>%
#   group_by(cluster) %>% summarise(cluster.region = unique(cluster.region))
#spring.gdata <-  merge(spring.gdata, spring.stat.regions, by = "cluster")

# Plot betweenness scores
ggplot(data = spring.gdata, aes(y = betweenness, x = as.factor(cluster), fill = node.type))+
  geom_col()+
  coord_flip()

# Plot bridge betweenness scores
ggplot(data = spring.gdata, aes(y = bridge.strength, x = as.factor(cluster), fill = node.type))+
  geom_col()+
  coord_flip()

## nonbreeding movement stats  ----

