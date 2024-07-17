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
library(jtools)
library(geosphere)
library(circular)

# Load the helper functions script  
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Load geolocator reference data 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# Load node data 
fall.ndata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/fall.graph.data.csv")
spring.ndata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/Spring.graph.data.csv")

#Load movement data
geo.all <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/All.locations.csv")

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

## migration distance ----

# Fall migration distance total
fall.dist <- fall.stat %>% group_by(geo_id) %>% mutate(Lon.50.next = lead(Lon.50.),
                                                       Lat.50.next = lead(Lat.50.)) %>%
  rowwise() %>%
  mutate(dists = distHaversine(c(Lon.50.,Lat.50.), c(Lon.50.next, Lat.50.next))) %>%
  ungroup() %>%
  group_by(geo_id) %>%
  filter(!is.na(dists)) %>%
  summarize(distance = sum(dists)/1000)

# Spring migration distance total 
spring.dist <- spring.stat %>% group_by(geo_id) %>% mutate(Lon.50.next = lead(Lon.50.),
                                                           Lat.50.next = lead(Lat.50.)) %>%
  rowwise() %>%
  mutate(dists = distHaversine(c(Lon.50.,Lat.50.), c(Lon.50.next, Lat.50.next))) %>%
  ungroup() %>%
  group_by(geo_id) %>%
  filter(!is.na(dists)) %>%
  summarize(distance = sum(dists)/1000)

# Fall migration between successive sites 
fall.dist.succ <- fall.stat %>% group_by(geo_id) %>% mutate(Lon.50.next = lead(Lon.50.),
                                                            Lat.50.next = lead(Lat.50.)) %>%
  rowwise() %>%
  mutate(dist.to.next.site = distHaversine(c(Lon.50.,Lat.50.), c(Lon.50.next, Lat.50.next)), .after = sitenum) 

# spring migration between successive sites 
spring.dist.succ <- spring.stat %>% group_by(geo_id) %>% mutate(Lon.50.next = lead(Lon.50.),
                                                                Lat.50.next = lead(Lat.50.)) %>%
  rowwise() %>%
  mutate(dist.to.next.site  = distHaversine(c(Lon.50.,Lat.50.), c(Lon.50.next, Lat.50.next)), .after = sitenum)

## Number of nodes used ----

# Load fall and spring inter node movement datasets
fall.edge.df.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.intra.cluster.movements.csv")
spring.edge.df.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.intra.cluster.movements.csv")

# load node data
meta.fall.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")
meta.spring.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")

# # remove supplementary edges in inter node movements (edges used in community analysis)
# mig.clusts.fall <- meta.fall.ab %>% dplyr::filter(node.type %in% c("Stopover", "Nonbreeding"))
# fall.edge.df.ab  <- fall.edge.df.ab %>% filter(next.cluster < max(mig.clusts.fall$vertex))
# 
# mig.clusts.spring <- meta.spring.ab %>% dplyr::filter(node.type %in% c("Stopover", "Nonbreeding"))
# spring.edge.df.ab  <- spring.edge.df.ab %>% filter(next.cluster < max(mig.clusts.spring$vertex))

# merge internode movements and node data to gert classification for each node used 
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

## function to calculate mode
Mode <- function(x) {
  ux <- unique(x) 
  ux <- ux[!(is.na(ux))]
  ux[which.max(tabulate(match(x, ux)))]
}

## number of nodes of each type 
meta.fall.ab  %>% group_by(node.type) %>% summarize(n())
meta.spring.ab  %>% group_by(node.type) %>% summarize(n())

### Average number of nodes used in fall network ----
range(analysis_ref$fall.nodes.occupied, na.rm =  T)
mean(analysis_ref$fall.nodes.occupied, na.rm =  T)
se <- sd(analysis_ref$fall.nodes.occupied, na.rm =  T)/sqrt(length(analysis_ref$fall.nodes.occupied[!is.na(analysis_ref$fall.nodes.occupied)]))
se
Mode(analysis_ref$fall.nodes.occupied)

### Average number of nodes used in the spring network  ----
range(analysis_ref$spring.nodes.occupied, na.rm =  T)
mean(analysis_ref$spring.nodes.occupied, na.rm =  T)
se <- sd(analysis_ref$spring.nodes.occupied, na.rm =  T)/sqrt(length(analysis_ref$spring.nodes.occupied[!is.na(analysis_ref$spring.nodes.occupied)]))
se
Mode(analysis_ref$spring.nodes.occupied)

### Average time spent stationary in each node ----
fall.node.times <- fall.stat %>% filter(!(is.na(duration))) %>% group_by(geo_id) %>%
  arrange(sitenum) %>%
  filter(sitenum != first(sitenum) & sitenum != last(sitenum)) %>% 
  group_by(geo_id, cluster) %>%
  summarize(mean.stat.dur = sum(duration)) 

mean(fall.node.times$mean.stat.dur)
se <- sd(fall.node.times$mean.stat.dur)/sqrt(length(fall.node.times$mean.stat.dur))
se
range(fall.node.times$mean.stat.dur)

spring.node.times <- spring.stat %>% filter(!(is.na(duration))) %>% group_by(geo_id) %>%
  arrange(sitenum) %>%
  filter(sitenum != first(sitenum) & sitenum != last(sitenum)) %>% ungroup() %>%
  group_by(geo_id, cluster) %>%
  summarize(mean.stat.dur = sum(duration)) 

mean(spring.node.times$mean.stat.dur)
se <- sd(spring.node.times$mean.stat.dur)/sqrt(length(spring.node.times$mean.stat.dur))
se
range(spring.node.times$mean.stat.dur)

### Duration of migration  ----

# Calculate the difference between the length of the fall and spring migration for every tracked individual 
migration_times <- geo.all %>% group_by(geo_id) %>% filter(StartTime == first(StartTime)) %>% summarize(fall.mig.duration = as.Date(nbr.arrival) - as.Date(fall.br.departure),
                                                                                                        spring.mig.duration = as.Date(spring.br.arrival) - as.Date(nbr.departure)) %>%
  mutate(timing.difference = fall.mig.duration - spring.mig.duration) %>%
  filter(!is.na(timing.difference))

# mean fall migration duration 
mean(migration_times$fall.mig.duration)
sd(migration_times$fall.mig.duration)/sqrt(length(migration_times$fall.mig.duration))

# mean spring migration duration 
mean(migration_times$spring.mig.duration)
sd(migration_times$spring.mig.duration)/sqrt(length(migration_times$spring.mig.duration))

# difference between spring and fall migration duration 
mean(migration_times$timing.difference)
sd(migration_times$timing.difference)/sqrt(length(migration_times$timing.difference))

#number of stopover nodes used by range region 
mean.stp.use.fall <- analysis_ref %>% group_by(Range_region) %>%
  summarize(mean.stp.used = mean(fall.stopover.nodes.occupied, na.rm = T), se.stp.used = mean.stp.used /sqrt(length(fall.stopover.nodes.occupied)),
            mode = Mode(fall.stopover.nodes.occupied))

mean.stp.use.spring <- analysis_ref %>% group_by(Range_region)%>%
  summarize(mean.stp.used = mean(spring.stopover.nodes.occupied, na.rm = T), se.stp.used = mean.stp.used /sqrt(length(spring.stopover.nodes.occupied)),
            mode = Mode(spring.stopover.nodes.occupied))

with(analysis_ref, table(fall.stopover.nodes.occupied, Range_region))
with(analysis_ref, table(spring.stopover.nodes.occupied, Range_region))

### correlation between number of  stopover nodes used and region of breeding site (eastern or western) ----
mod.fall.stopover.node <- glmmTMB(fall.stopover.nodes.occupied   ~ as.factor(Range_region), data = analysis_ref, family = genpois(link = "log"))
boxplot(fall.stopover.nodes.occupied ~ Range_region, data = analysis_ref)
summary(mod.fall.stopover.node)
#check_model(mod.fall.stopover.node)

simulationOutput <- simulateResiduals(fittedModel = mod.fall.stopover.node, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

mod.spring.stopover.node <- glmmTMB(spring.stopover.nodes.occupied  ~ as.factor(Range_region), data = analysis_ref, family = genpois(link = "log"))
boxplot(spring.stopover.nodes.occupied ~ Range_region, data = analysis_ref)
summary(mod.spring.stopover.node)
#check_model(mod.spring.stopover.node)

simulationOutput <- simulateResiduals(fittedModel = mod.spring.stopover.node , plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

## Number of nonbreeding node used ----

#number of nonbreeding nodes used by range region 
mean.nbr.use.fall <- analysis_ref %>% group_by(Range_region) %>%
  summarize(mean.nbr.used = mean(fall.nbr.nodes.occupied, na.rm = T), se.nbr.used = mean.nbr.used /sqrt(length(fall.nbr.nodes.occupied)),
            mode = Mode(fall.nbr.nodes.occupied))

mean.nbr.use.spring <- analysis_ref %>% group_by(Range_region)%>%
  summarize(mean.nbr.used = mean(spring.nbr.nodes.occupied, na.rm = T), se.nbr.used = mean.nbr.used /sqrt(length(spring.nbr.nodes.occupied)),
            mode = Mode(spring.nbr.nodes.occupied))

with(analysis_ref, table(fall.nbr.nodes.occupied, Range_region))
with(analysis_ref, table(spring.nbr.nodes.occupied, Range_region))

### correlation between number of  nonbreeding nodes used and region of breeding site (eastern or western) ----
mod.fall.nbr.node <- glmmTMB(fall.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref, family = genpois(link = "log"))
boxplot(fall.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref)
summary(mod.fall.nbr.node)
#check_model(mod.fall.nbr.node)

simulationOutput <- simulateResiduals(fittedModel = mod.fall.nbr.node, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

mod.spring.nbr.node <- glm(spring.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref, family = poisson(link = "log"))
boxplot(spring.nbr.nodes.occupied ~ as.factor(Range_region), data = analysis_ref)
summary(mod.spring.nbr.node)
#check_model(mod.spring.nbr.node)

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
  scale_colour_manual(values = c("Northwestern Region" = "#F0E442",
                               "Central Region" = "#009E73",
                               "Eastern Region" = "#0072B2",
                               "Western Region"  = "#D55E00"), name = "Breeding region") +
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "First nonbreeding site longitude") # +
#geom_text(label = geo.loc.data$geo_id, nudge_x = 5, nudge_y = 1)

# run a regression model
mod1.data <- analysis_ref %>% filter_at(vars(deploy.longitude, deploy.latitude, nbr1.lon),all_vars(!is.na(.)))
mod1 <- glmmTMB(nbr1.lon ~ deploy.longitude +(1|study.site), data = mod1.data, na.action = "na.fail")
# # mod1 <- plsr(nbr1.lon ~ deploy.latitude + deploy.longitude, data = mod1.data, scale = T, validation = "CV",
# #              ncomp = 1, method = "oscorespls")
# mod1 <- glmmTMB(nbr1.lon ~ poly(deploy.longitude, degree = 2) + (1|study.site), data = mod1.data, na.action = "na.fail")
summary(mod1)
# check_model(mod1)
# #with(mod1.data, table(study.site))
# 
simulationOutput <- simulateResiduals(fittedModel =  mod1, plot = F, quantreg = T)
plot(simulationOutput)
plotResiduals(simulationOutput)
# 
effect_plot(mod1, pred = deploy.longitude, interval = TRUE, plot.points = TRUE, 
            x.label = "Breeding longitude",
            y.label = "nonbreeding longitude")

# spearman's rank correlation
cor.test(mod1.data$nbr1.lon, mod1.data$deploy.longitude, method = "spearman", exact = F)

## Plot the second nonbreeding longitude against the breeding longitude ----
long.plot2 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr2.lon)) + 
  geom_point(aes(col = Breeding_region_MC)) +
  scale_colour_manual(values = c("Northwestern Region" = "#F0E442",
                                 "Central Region" = "#009E73",
                                 "Eastern Region" = "#0072B2",
                                 "Western Region"  = "#D55E00"), name = "Breeding region") +
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme(text = element_text(size=14)) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding longitude",
       y = "second nonbreeding site longitude")

# run a regression model
mod2.data <- analysis_ref %>% filter_at(vars(deploy.longitude, deploy.latitude, nbr2.lon),all_vars(!is.na(.)))
mod2 <- glmmTMB(nbr2.lon ~ deploy.longitude + (1|study.site), data = mod2.data, na.action = "na.fail")
# mod2 <- lm(nbr2.lon ~ deploy.longitude, data = mod2.data, na.action = "na.fail")
summary(mod2)
# check_model(mod2)
# 
# # mod1.data.reg1 <- mod1.data %>% filter(Breeding_region_MC == "Eastern Region")
# # cor(mod1.data.reg1$deploy.latitude, mod1.data.reg1$deploy.longitude)
# # mod1.data.reg2 <- mod1.data %>% filter(Breeding_region_MC == "Central Region")
# # plot(mod1.data.reg2$deploy.latitude, mod1.data.reg2$deploy.longitude)
# # mod1.data.reg3 <- mod1.data %>% filter(Breeding_region_MC == "Western Region")
# # plot(mod1.data.reg3$deploy.latitude, mod1.data.reg3$deploy.longitude)
# # mod1.data.reg4 <- mod1.data %>% filter(Breeding_region_MC == "Northwestern Region")
# # cor(mod1.data.reg4$deploy.latitude, mod1.data.reg4$deploy.longitude)
# 
simulationOutput <- simulateResiduals(fittedModel =  mod2, plot = F)
plot(simulationOutput)

effect_plot(mod2, pred = deploy.longitude, interval = TRUE, plot.points = TRUE, 
            x.label = "Breeding longitude",
            y.label = "nonbreeding longitude")

# spearman's rank correlation
cor.test(mod2.data$nbr2.lon, mod2.data$deploy.longitude, method = "spearman", exact = F)

## Plot the first nonbreeding latitude against the breeding latitude ----
lat.plot1 <- ggplot(data = analysis_ref, aes(x = deploy.latitude, y = nbr1.lat)) + 
  geom_point(aes(col = Breeding_region_MC)) + 
  scale_colour_manual(values = c("Northwestern Region" = "#F0E442",
                                 "Central Region" = "#009E73",
                                 "Eastern Region" = "#0072B2",
                                 "Western Region"  = "#D55E00"), name = "Breeding region") +
  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "First nonbreeding site latitude")

# run a regression model
mod3.data <- analysis_ref %>% filter_at(vars(deploy.latitude, nbr1.lat),all_vars(!is.na(.)))
mod3 <- glmmTMB(nbr1.lat ~ deploy.latitude + (1|study.site), data = mod3.data )
# plot(nbr1.lat ~ deploy.latitude, data = mod3.data )
summary(mod3)
# check_model(mod3)
# 
simulationOutput <- simulateResiduals(fittedModel =  mod3, plot = F)
plot(simulationOutput)

effect_plot(mod3, pred = deploy.latitude, interval = TRUE, plot.points = TRUE, 
            x.label = "Breeding latitude",
            y.label = "Nonbreeding latitude")
            
# 
## Plot the second nonbreeding latitude against the breeding latitude ----
lat.plot2 <- ggplot(data = analysis_ref, aes(x = deploy.latitude, y = nbr2.lat)) +
  geom_point(aes(col = Breeding_region_MC)) +
  scale_colour_manual(values = c("Northwestern Region" = "#F0E442",
                                 "Central Region" = "#009E73",
                                 "Eastern Region" = "#0072B2",
                                 "Western Region"  = "#D55E00"), name = "Breeding region") +

  geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Breeding latitude",
       y = "Second nonbreeding site latitude")

# spearman's rank correlation
cor.test(mod3.data$nbr1.lat, mod3.data$deploy.latitude, method = "spearman", exact = F)

# run a regression model
mod4.data <- analysis_ref %>% filter_at(vars(deploy.latitude, nbr2.lat),all_vars(!is.na(.)))
mod4 <- glmmTMB(nbr2.lat ~ deploy.latitude + (1|study.site), data = mod4.data, na.action = "na.fail")
# mod4 <- lm(nbr2.lat ~ deploy.latitude, data = mod4.data )
# plot(nbr2.lat ~ deploy.latitude, data = mod4.data )
summary(mod4)
# check_model(mod4)
# 
simulationOutput <- simulateResiduals(fittedModel =  mod4, plot = F)
plot(simulationOutput)

effect_plot(mod4, pred = deploy.latitude, interval = TRUE, plot.points = TRUE, 
            x.label = "Breeding latitude",
            y.label = "Nonbreeding latitude")
# 
# x <- analysis_ref$nbr2.lon[!is.na(analysis_ref$nbr2.lon)]
# y <- analysis_ref$deploy.longitude[!is.na(analysis_ref$nbr2.lon)]
# cor(x,y, method = "pearson")
# 
# ## Plot the first nonbreeding latitude against the breeding longitude ----
# lonlat.plot1 <- ggplot(data = analysis_ref, aes(x = deploy.longitude, y = nbr1.lat)) + 
#   geom_point() + 
#   geom_smooth(method = "lm", se = F, colour="black", size=0.5, linetype = "dashed") + 
#   theme_bw() +
#   theme(text = element_text(size = 14)) +
#   labs(x = "Breeding longitude",
#        y = "First nonbreeding site latitude")

# spearman's rank correlation
cor.test(mod4.data$nbr2.lat, mod4.data$deploy.latitude, method = "spearman", exact = F)

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
NB.stat.mean <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Data/nonbreeding.mean.csv")

# Update reference and location data with the status of each individual: whether or not they performed nonbreeding movements 
analysis_ref <- analysis_ref %>% mutate(move.status = ifelse(geo_id %in% NB.move$geo_id, "mover", "nonmover")) %>%
  filter(geo_id %in% NB.move$geo_id | geo_id %in% NB.stat.mean$geo_id)
geo.all.nbr <- geo.all %>% mutate(move.status = ifelse(geo_id %in% NB.move$geo_id, "mover", "nonmover")) %>%
  filter(geo_id %in% NB.move$geo_id | geo_id %in% NB.stat.mean$geo_id)

## How many birds performed nonbreeding movements? ----
analysis_ref %>% group_by(move.status) %>% summarize(n = n())
#View(NB.mover.cat %>% group_by(status) %>% summarize(geo_id = unique(geo_id)))

## test whether the probability of nonbreeding movements was influenced by the longitude of the first nonbreeding site occupied ---- 
nbr.locations <- geo.all.nbr  %>% filter(NB_count == 1)  
nbr.lon.mod <- glm(as.factor(move.status) ~ Lon.50., data = nbr.locations, family = binomial(link = "logit"))
boxplot(as.factor(move.status) ~ Lon.50., data = nbr.locations)
summary(nbr.lon.mod)
check_model(nbr.lon.mod)

simulationOutput <- simulateResiduals(fittedModel = nbr.lon.mod, plot = F)
plot(simulationOutput)
testDispersion(simulationOutput)

# movement timing summary 
NB.move %>% group_by(timing.nbr.move) %>%summarise (n = n())

# movement occurrence by month 
NB.move <- NB.move %>% mutate(move.month = month(move.start))
NB.move %>% group_by(move.month) %>% summarise(n())

# movement direction (North south) summary
NB.move <- NB.move %>% rowwise() %>% mutate(move.bearing = bearing(c(start.lon, start.lat), c(end.lon, end.lat)), 
                                            move.direction = ifelse(abs(move.bearing) < 90, "north", "south")) 
print(NB.move %>% group_by(move.direction) %>% summarise (n = length(unique(geo_id))), n = 24)

# time spent at northern nonbreeding grounds 
northward.nbr <- NB.move %>% group_by(geo_id) %>% filter(move.direction == "north", StartTime == max(StartTime),
                                                         duration < 100)
mean(northward.nbr$duration) 
range(northward.nbr$duration)
sd(northward.nbr$duration)/sqrt(length(northward.nbr$duration))

# time spent at southern nonbreeding grounds (and one nonbreeding grounds from a northward movement in the early winter)
southward.nbr <- NB.move %>% group_by(geo_id) %>% filter(move.direction == "south" | duration > 100)

mean(southward.nbr$duration) 
range(southward.nbr$duration)
sd(southward.nbr$duration)/sqrt(length(southward.nbr$duration))

# individual bird movement directions
move.schedule <- NB.move %>% group_by(geo_id,move.direction) %>%
  reframe(n = n(), move.start) %>%
  arrange(move.start)
print(move.schedule, n = nrow(move.schedule))

# movement proximity to the equinox summamove.month
NB.move %>% 
  group_by(move.direction, equinox.nbr.move) %>%summarise (n = n())

#time spent at the last nonbreeding sites occupied by individuals that moved
geo.all %>% group_by(geo_id) %>% filter(geo_id %in% NB.move$geo_id, NB_count == max(NB_count, na.rm = T)) %>% 
  summarize(time.last.site = max(duration)) 

# Average nonbreeding movement distance mean and SE ----
mean(NB.move$dist)
sd(NB.move$dist)/sqrt(length(NB.move$dist))

## Statistics for the number of movements per bird ----
NB.move %>%  group_by(geo_id)%>% summarize(move.num = n()) %>% group_by(move.num) %>%
  summarize(birds = n())

## Statistics for the bearing of movements between the first and last nonbreeding sites ---- 
NB.all <- geo.all %>% group_by(geo_id) %>% filter(NB_count == min(NB_count, na.rm =  T) | NB_count == max(NB_count, na.rm =  T),
                                                  geo_id %in% NB.move$geo_id) %>%
  mutate(Lon.50.next = lead(Lon.50.),
         Lat.50.next = lead(Lat.50.)) %>% rowwise() %>%
  mutate(bearing = bearing(c(Lon.50., Lat.50.), c(Lon.50.next, Lat.50.next)),
         course = (bearing + 360) %% 360) %>%
  filter(!is.na(bearing)) %>% dplyr::select(geo_id, bearing, course)

mean(circular(NB.all$course, units = "degrees"))
rayleigh.test(circular(NB.all$bearing, units = "degrees"))

nrow(NB.all[abs(NB.all$bearing) < 90,])
nrow(NB.all)

# ggplot(st_as_sf(wrld_simpl))+
#   geom_sf(colour = "black", fill = "#F7F7F7", lwd = 0.3) +
#   coord_sf(xlim = c(-95, -45),ylim = c(-10, 15)) +
#   geom_arrowsegment(data = NB.all[month(NB.all$EndTime) %in% c(2,3,4, 5),], mapping = aes(x = Lon.50., y = Lat.50., xend = Lon.50.next, yend = Lat.50.next),
#                     arrows = arrow(end = "last", type = "closed", length = unit(0.1, "inches")), arrow_positions = 1, lwd = 0.6)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
#         plot.title=element_text(size=8, vjust=-1),
#         legend.position = c(0.14, 0.4),
#         axis.title =element_blank(),
#         axis.text =element_blank(),
#         axis.ticks =element_blank(),
#         axis.ticks.length = unit(0, "pt"),
#         legend.spacing = unit(-5, "pt"),
#         plot.margin = unit(c(0,0,0,0), "pt"),
#         legend.key = element_rect(colour = "transparent", fill = "white"))

# Average arrival and departure dates in the nonbreeding grounds ----

# average date of arrival from the nonbreeding grounds 
avg.nbr.arr <- mean(yday(geo.all$nbr.arrival), na.rm = T)
as.Date(avg.nbr.arr, origin = "2020-01-01")

# average date of departure from the nonbreeding grounds 
avg.nbr.dep <- mean(yday(geo.all$nbr.departure), na.rm = T)
as.Date(avg.nbr.dep, origin = "2020-01-01")

# Transoceanic crossings ----

# shapefile with fall stationary locations
fall.stat.sf <- st_as_sf(fall.stat, coords = c("Lon.50.", "Lat.50."), crs = crs(wrld_simpl))

# load Caribbean polygon 
caribb.poly <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Mapping_components/Data/Caribb.poly.shp") %>%
  st_transform(st_crs(fall.stat.sf))

caribb.countries <- st_as_sf(wrld_simpl[wrld_simpl$SUBREGION == 29,])

# get the fall caribbean stopovers 
fall.carib.stops1 <- st_intersection(fall.stat.sf , caribb.countries)
fall.carib.stops2 <- st_intersection(fall.stat.sf , caribb.poly)

# convert fall locations to spatial points 
fall.stat.sf <- st_as_sf(fall.stat, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl), remove = F)

# convert spring locations to spatial points 
spring.stat.sf <- st_as_sf(spring.stat, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl), remove = F)

# get the spring caribbean stopovers 
spring.carib.stops1 <- st_intersection(spring.stat.sf , caribb.countries)
spring.carib.stops2 <- st_intersection(spring.stat.sf , caribb.poly)

# plot of carribean stopovers
America <- wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),]
ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = fall.carib.stops2, colour = "#D55E00")+
  geom_sf(data = spring.carib.stops2, colour = "#009E73")+
  #geom_sf_text(data = fall.carib.stops2, aes(label =  geo_id))+
  #geom_sf_text(data = spring.carib.stops2, aes(label =  geo_id))+
  coord_sf(xlim = c(-90, -50),ylim = c(10, 30)) 


# Add stopover information to reference data 
analysis_ref_carib <- analysis_ref %>% mutate(fall.carib.stop = ifelse(geo.id %in% fall.carib.stops2$geo_id, "stopover", "direct"),
                                              spring.carib.stop = ifelse(geo.id %in% spring.carib.stops2$geo_id, "stopover", "direct"),
                                              fall.carib.stop = ifelse(geo.id %in% fall.stat$geo_id, fall.carib.stop, NA),
                                              spring.carib.stop = ifelse(geo.id %in% spring.stat$geo_id, spring.carib.stop, NA))

x <- analysis_ref_carib %>% group_by(Breeding_region_MC, fall.carib.stop) %>% summarize(unique(geo.id)) 


# Estimate geolocator error and bias ----

# the projection used for location data 
proj <- '+proj=aeqd +lat_0=0 +lon_0=-74'


# Get locations obtained by runnign the threshold group model without grouping, projected
locs <- findThresModData() 
locs <- locs %>% group_by(geo_id) %>% filter(Time1 > IH.calib.start & Time1  < IH.calib.end)


Thresh.mod.data.loc <- locs %>% st_as_sf(coords = c("Lon.mean",
                                                    "Lat.mean"),
                                         crs = 4326) %>%
  st_transform(proj)

# origin locations, projected
Thresh.mod.data.or <- locs %>% st_as_sf(coords = c("deploy.longitude",
                                                   "deploy.latitude"),
                                        crs = 4326) %>%
  st_transform(proj)

# geolocator points 
geo.br <- rbind(fall.breed, spring.breed[!(spring.breed$geo_id %in% fall.breed$geo_id),]) %>% arrange(geo_id)
geo.br.sf <- st_as_sf(geo.br, coords = c("deploy.longitude", "deploy.latitude"), crs = crs(wrld_simpl))
geo.br.sf <- st_transform(geo.br.sf, CRS(proj)) 

# We calculate spatial bias as the mean distance from the geolocator deployment site while the bird was known to be at that location
lon.errors <- c()
lat.errors <- c()
lon.biass <- rep(NA, nrow(geo.br.sf ))
lat.biass <- rep(NA, nrow(geo.br.sf ))
lon.sds <- rep(NA, nrow(geo.br.sf ))
lat.sds <- rep(NA, nrow(geo.br.sf ))

for (i in seq(1, length(geo.br.sf$geo_id))){
  
  mod.locs <- Thresh.mod.data.loc %>% filter(geo_id == geo.br.sf$geo_id[i])
  or.locs <-  Thresh.mod.data.or %>% filter(geo_id == geo.br.sf$geo_id[i])
  
  lon.bias <- mean(st_coordinates(mod.locs)[,1] - st_coordinates(or.locs)[,1])
  lat.bias <- mean(st_coordinates(mod.locs)[,2] - st_coordinates(or.locs)[,2])
  lon.error <- abs(st_coordinates(mod.locs)[,1] - st_coordinates(or.locs)[,1])
  lat.error <- abs(st_coordinates(mod.locs)[,2] - st_coordinates(or.locs)[,2])
  
  lon.sd <- sd(st_coordinates(mod.locs)[,1])
  lat.sd <- sd(st_coordinates(mod.locs)[,2])
  
  lon.biass[i] <- lon.bias
  lat.biass[i] <- lat.bias
  lon.errors <- append(lon.errors, lon.error)
  lat.errors <- append(lat.errors, lat.error)
  lon.sds[i] <- lon.sd
  lat.sds[i] <- lat.sd
  lon.sds[i] <- lon.sd
  lat.sds[i] <- lat.sd
  
}

# Geolocator error and variance/covariance at origin sites
mod <- lm(cbind(lon.biass, lat.biass) ~ 1)
geo.bias <- coef(mod)
geo.vcov <- vcov(mod)

mean(lon.errors)
sd(lon.errors)
mean(lat.errors)
sd(lat.errors)

# Assesss the proportion of individuals using a node 

spring.stat %>% filter(geo_id %in% spring.stat[spring.stat$cluster %in% c(11),]$geo_id) %>%
  group_by(Breeding_region_MC) %>%
  summarize(n = length(unique(geo_id)))

spring.stat %>% group_by(Breeding_region_MC) %>%
  summarize(n = length(unique(geo_id)))

spring.stat %>% filter(!(geo_id %in% spring.stat[(spring.stat$cluster %in% c(11)),]$geo_id)) %>%
  group_by(Breeding_region_MC) %>%
  summarize(n = length(unique(geo_id)))

spring.stat.ab %>% filter((geo_id %in% spring.stat[(spring.stat$cluster %in% c(11)),]$geo_id)) %>%
  group_by(Breeding_region_MC, geo_id) %>%
  summarize(ab.unit = unique(ab.unit)) %>%
  summarize(ab.total = sum(ab.unit))

fall.stat %>% filter(geo_id %in% fall.stat[fall.stat$cluster %in% c(3,5),]$geo_id) %>%
  group_by(Breeding_region_MC) %>%
  summarize(n = length(unique(geo_id)))

fall.stat.ab %>% filter(geo_id %in% fall.stat.ab[fall.stat.ab$cluster %in% c(5),]$geo_id) %>%
  group_by(Breeding_region_MC, geo_id) %>%
  summarize(ab.unit = unique(ab.unit)) %>%
  summarize(ab.total = sum(ab.unit))