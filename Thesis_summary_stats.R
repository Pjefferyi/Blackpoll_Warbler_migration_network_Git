# Thesis summary stats  

# load necessary libraries 
library(performance)
library(DHARMa)
library(MASS)
library(tidyr)
library(purrr)
library(glmmTMB)
library(MuMIn)
library(lme4)

# Load tMuMIn# Load the helper functions script  
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Load geolocator reference data 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# Load node data 
fall.ndata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/fall.graph.data.csv")
spring.ndata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/Spring.graph.data.csv")

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
              "E",
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
tol_values <- seq(0, 0.26, by = 0.001) 
fall_equi <- anytime("2023-09-22")
cens_data <- data.frame(tol_values)

cens_data <- cens_data %>% 
  rowwise() %>%
  mutate(days_cens = length(solar.angle[solar.angle > 0 - tol_values & solar.angle < 0 + tol_values])/2)

# calculate the number of days censored for each geolocator 
analysis_ref <- analysis_ref %>% rowwise() %>% mutate(tol_days = list(cens_data[round(tol_values, digits = 3) == round(tol, digits = 3),]$days_cens))
analysis_ref$tol_days <- unlist(analysis_ref$tol_days)

# Range of days 
range(analysis_ref$tol_days)

# Mean and SD
mean(analysis_ref$tol_days)
sd(analysis_ref$tol_days)

# Results #####################################################################

# Load fall and spring stationary sites 
fall.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.data.csv") %>%
  merge(fall.ndata, by.x = "cluster", by.y ="X") 
spring.stat <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.stationary.data.csv") %>%
  merge(fall.ndata, by.x = "cluster", by.y ="X") 

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
fall.node.used.df <- fall.edge.df.ab %>% group_by(geo_id) %>% reframe(fall.nodes.occupied = length(unique(cluster, next.cluster)),
                                                                        fall.stopover.nodes.occupied = length(unique(cluster[node.type == "Stopover"], next.cluster[node.type == "Stopover"])),
                                                                        fall.nbr.nodes.occupied = length(unique(cluster[node.type == "Nonbreeding"], next.cluster[node.type == "Nonbreeding"]))) 

# number of nodes used in the spring network 
spring.node.used.df <- spring.edge.df.ab %>% group_by(geo_id) %>% reframe(spring.nodes.occupied = length(unique(cluster, next.cluster)),
                                                                        spring.stopover.nodes.occupied = length(unique(cluster[node.type == "Stopover"], next.cluster[node.type == "Stopover"])),
                                                                        spring.nbr.nodes.occupied = length(unique(cluster[node.type == "Nonbreeding"], next.cluster[node.type == "Nonbreeding"])))

# merge the information collected with the reference dataset 
analysis_ref <- merge(analysis_ref, fall.node.used.df, by.x = "geo.id", by.y = "geo_id", all = T)
analysis_ref <- merge(analysis_ref, spring.node.used.df, by.x = "geo.id", by.y = "geo_id", all = T) %>%
  mutate(Range_region = ifelse(Range_region == "Central", "Western", Range_region))

## Statistics of node usage ----

## number of nodes of each type 
meta.fall.ab  %>% group_by(node.type) %>% summarize(n())
meta.spring.ab  %>% group_by(node.type) %>% summarize(n())

### Average number of nodes used in fall network ----
range(analysis_ref$fall.nodes.occupied, na.rm =  T)
mean(analysis_ref$fall.nodes.occupied, na.rm =  T)
se <- sd(analysis_ref$fall.nodes.occupied, na.rm =  T)/sqrt(length(analysis_ref$fall.nodes.occupied[!is.na(analysis_ref$fall.nodes.occupied)]))
se

### Average number of nodes used in the spring network  ----
range(analysis_ref$spring.nodes.occupied, na.rm =  T)
mean(analysis_ref$spring.nodes.occupied, na.rm =  T)
se <- sd(analysis_ref$spring.nodes.occupied, na.rm =  T)/sqrt(length(analysis_ref$spring.nodes.occupied[!is.na(analysis_ref$spring.nodes.occupied)]))
se

# ### correlation between number of  nodes used and longitude of breeding site ----
# mod.fall.stopover.node <- glm(fall.stopover.nodes.occupied   ~ deploy.longitude, data = analysis_ref, family = gaussian(link = "identity"))
# plot(fall.stopover.nodes.occupied ~ deploy.longitude, data = analysis_ref)
# summary(mod.fall.stopover.node)
# check_model(mod.fall.stopover.node)
# 
# simulationOutput <- simulateResiduals(fittedModel = mod.fall.stopover.node, plot = F)
# plot(simulationOutput)
# 
# mod.spring.stopover.node <- glm(spring.stopover.nodes.occupied  ~ (deploy.longitude), data = analysis_ref, family = gaussian(link = "identity"))
# plot(spring.stopover.nodes.occupied  ~ deploy.longitude, data = analysis_ref)
# summary(mod.spring.stopover.node)
# check_model(mod.spring.stopover.node)
# 
# simulationOutput <- simulateResiduals(fittedModel = mod.spring.stopover.node , plot = F)
# plot(simulationOutput)

### correlation between number of  stopovernodes used and region of breeding site (eastern or western) ----
mod.fall.stopover.node <- glmmTMB(fall.stopover.nodes.occupied   ~ as.factor(Range_region), data = analysis_ref, family = genpois(link = "log"))
boxplot(fall.stopover.nodes.occupied ~ Range_region, data = analysis_ref)
summary(mod.fall.stopover.node)
check_model(mod.fall.stopover.node)

simulationOutput <- simulateResiduals(fittedModel = mod.fall.stopover.node, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

mod.spring.stopover.node <- glmmTMB(spring.stopover.nodes.occupied  ~ as.factor(Range_region), data = analysis_ref, family = genpois(link = "log"))
boxplot(spring.stopover.nodes.occupied ~ Range_region, data = analysis_ref)
summary(mod.spring.stopover.node)
check_model(mod.spring.stopover.node)

simulationOutput <- simulateResiduals(fittedModel = mod.spring.stopover.node , plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

### correlation between number of  stopovernodes used and region of breeding site (eastern or western) ----
mod.fall.nbr.node <- glmmTMB(fall.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref, family = genpois(link = "log"))
boxplot(fall.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref)
summary(mod.fall.nbr.node)
check_model(mod.fall.nbr.node)

simulationOutput <- simulateResiduals(fittedModel = mod.fall.nbr.node, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

mod.spring.nbr.node <- glm(spring.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref, family = poisson(link = "log"))
boxplot(spring.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref)
summary(mod.spring.nbr.node)
check_model(mod.spring.nbr.node)

simulationOutput <- simulateResiduals(fittedModel = mod.spring.stopover.node, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

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
analysis_ref <- loc_list %>% purrr::reduce(full_join, by='geo_id') %>% 
  merge(analysis_ref, by.x = "geo_id", by.y = "geo.id", all = T)

## Plot the first nonbreeding longitude against the breeding longitude ----
long.plot1 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr1.lon)) + 
  geom_point(aes(col = Breeding_region_MC)) +
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "First nonbreeding site longitude") # +
#geom_text(label = geo.loc.data$geo_id, nudge_x = 5, nudge_y = 1)

# run a regression model
mod1.data <- analysis_ref %>% filter_at(vars(deploy.longitude, deploy.latitude, nbr1.lon),all_vars(!is.na(.)))
#mod1 <- glmmTMB(nbr1.lon ~ deploy.latitude + (1|study.site), data = mod1.data, na.action = "na.fail")
# mod1 <- plsr(nbr1.lon ~ deploy.latitude + deploy.longitude, data = mod1.data, scale = T, validation = "CV",
#              ncomp = 1, method = "oscorespls")
mod1 <- lm(nbr1.lon ~ deploy.longitude, data = mod1.data, na.action = "na.fail")
summary(mod1)
check_model(mod1)
#with(mod1.data, table(study.site))

simulationOutput <- simulateResiduals(fittedModel =  mod1, plot = F, quantreg = T)
plot(simulationOutput)

## Plot the second nonbreeding longitude against the breeding longitude ----
long.plot2 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr2.lon)) + 
  geom_point(aes(col = Breeding_region_MC)) +
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme(text = element_text(size=14)) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "second nonbreeding site longitude")

# run a regression model
mod2.data <- analysis_ref %>% filter_at(vars(deploy.longitude, deploy.latitude, nbr2.lon),all_vars(!is.na(.)))
#mod2 <- glmmTMB(nbr2.lon ~ deploy.longitude + deploy.latitude + (1|study,site), data = mod2.data, na.action = "na.fail")
mod2 <- lm(nbr2.lon ~ deploy.longitude, data = mod2.data, na.action = "na.fail")
summary(mod2)
check_model(mod2)

# mod1.data.reg1 <- mod1.data %>% filter(Breeding_region_MC == "Eastern Region")
# cor(mod1.data.reg1$deploy.latitude, mod1.data.reg1$deploy.longitude)
# mod1.data.reg2 <- mod1.data %>% filter(Breeding_region_MC == "Central Region")
# plot(mod1.data.reg2$deploy.latitude, mod1.data.reg2$deploy.longitude)
# mod1.data.reg3 <- mod1.data %>% filter(Breeding_region_MC == "Western Region")
# plot(mod1.data.reg3$deploy.latitude, mod1.data.reg3$deploy.longitude)
# mod1.data.reg4 <- mod1.data %>% filter(Breeding_region_MC == "Northwestern Region")
# cor(mod1.data.reg4$deploy.latitude, mod1.data.reg4$deploy.longitude)

simulationOutput <- simulateResiduals(fittedModel =  mod2, plot = F)
plot(simulationOutput)

## Plot the first nonbreeding latitude against the breeding latitude ----
lat.plot1 <- ggplot(data = analysis_ref, aes(x = deploy.latitude, y = nbr1.lat)) + 
  geom_point(aes(col = Breeding_region_MC)) + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "First nonbreeding site latitude")

# run a regression model
mod3.data <- analysis_ref %>% filter_at(vars(deploy.latitude, nbr1.lat),all_vars(!is.na(.)))
mod3 <- lm(nbr1.lat ~ deploy.longitude, data = mod3.data )
plot(nbr1.lat ~ deploy.longitude, data = mod3.data )
summary(mod3)
check_model(mod3)

simulationOutput <- simulateResiduals(fittedModel =  mod3, plot = F)
plot(simulationOutput)

## Plot the second nonbreeding latitude against the breeding latitude ----
lat.plot2 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr2.lat)) + 
  geom_point(aes(col = Breeding_region_MC)) + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "Second nonbreeding site latitude")

# run a regression model
mod4.data <- analysis_ref %>% filter_at(vars(deploy.latitude, nbr2.lat),all_vars(!is.na(.)))
#mod4 <- glmmTMB(nbr2.lat ~ deploy.latitude + (1|study_site), data = mod4.data, na.action = "na.fail")
mod4 <- lm(nbr2.lat ~ deploy.longitude, data = mod4.data )
plot(nbr2.lat ~ deploy.longitude, data = mod4.data )
summary(mod4)
check_model(mod4)

simulationOutput <- simulateResiduals(fittedModel =  mod4, plot = F)
plot(simulationOutput)

x <- analysis_ref$nbr2.lon[!is.na(analysis_ref$nbr2.lon)]
y <- analysis_ref$deploy.longitude[!is.na(analysis_ref$nbr2.lon)]
cor(x,y, method = "pearson")

## Plot the first nonbreeding latitude against the breeding longitude ----
lonlat.plot1 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr1.lat)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "First nonbreeding site latitude")

## network metric scores by nodes ----

#load network data
fall.gdata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/fall.graph.data.csv")
spring.gdata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/spring.graph.data.csv")

# merge with fall graph data with fall.stat region names 
fall.stat.clusters <- fall.stat %>%
  group_by(cluster) %>% summarise(cluster = unique(cluster))
fall.gdata <-  merge(fall.gdata, fall.stat.clusters, by = "cluster")

# Plot betweenness scores
ggplot(data = fall.gdata, aes(y = betweenness, x = factor(cluster), fill = node.type))+
  geom_col()+
  coord_flip()
  
# Plot indegree strengthscores
ggplot(data = fall.gdata, aes(y = bridge.indegree, x = factor(cluster), fill = node.type))+
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

# Plot bridge indegree strength scores
ggplot(data = spring.gdata, aes(y = bridge.indegree, x = as.factor(cluster), fill = node.type))+
  geom_col()+
  coord_flip()

## nonbreeding movement stats  ----

# read nonbreeding movement data generated with this script C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Blackpoll_warbler_mapping_scripts/Blackpoll_nonbreeding_movements.R
NB.move <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/nonbreeding.movements.csv")

# get dataframe with categorical information on movement vs no-movement, longituide of the first nonbreeding site, and number of movements per individual
NB.mover.cat <- NB.move %>% group_by(geo_id) %>%
  mutate(initial.nbr.lon = first(Lon.50.))%>%
  summarise(status = unique(nbr.mover), n.movements = n(), initial.nbr.lon = first(initial.nbr.lon), breeding.lon = first(deploy.longitude)) %>%
  mutate(n.movements = ifelse(status == "nonmover", 0, n.movements/2))
  
## How many birds performed nonbreeding movements? ----
NB.mover.cat %>% group_by(status) %>% summarize(n = n())

## test whether the probability of nonbreeding movements was influenced by the longitude of the first nonbreeding site occupied ---- 
nbr.lon.mod <- glm(as.factor(status) ~ initial.nbr.lon, data = NB.mover.cat, family = binomial(link = "logit"))
boxplot(initial.nbr.lon ~ as.factor(status), data = NB.mover.cat)
summary(nbr.lon.mod)
check_model(nbr.lon.mod)

## get dataframe for movements that classifies timing and direction -----
NB.move.mod <- NB.move %>% filter(nbr.mover == "mover") %>% group_by(geo_id) %>%
  mutate(move.start = EndTime,
         move.end = lead(StartTime),
         start.lon = Lon.50.,
         start.lat = Lat.50.,
         end.lon = lead(Lon.50.),
         end.lat = lead(Lat.50.),
         dist = lead(dist)) %>% dplyr::select(geo_id, move.start, move.end, start.lon, start.lat, end.lon, end.lat, dist, timing.nbr.move)%>%
  filter(move.start < move.end, !is.na(move.end)) %>%
  mutate(move.direction = ifelse(start.lat < end.lat, "North", "South"))

NB.move.mod %>% group_by(timing.nbr.move) %>%summarise (n = n())
NB.move.mod %>% group_by(move.direction) %>%summarise (n = n())


