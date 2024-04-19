# clear objects from workspace
#rm(list=ls())

#Load libraries
library(tidyverse)
library(igraph)
library(lubridate)
library(anytime)
library(ggplot2)
library(geosphere)
library(terra)
library(sf)
library(maptools)
library(ebirdst)
library(scatterpie)

source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Import data ##################################################################

# path to reference data file 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# location data 
geo.all <- findLocData(geo.ids = c("V8757_010",
                                   "V8296_004",
                                   "V8296_005",
                                   "V8296_006",
                                   "V8757_055",
                                   "V8757_018",
                                   "V8296_015",
                                   "V8296_017",
                                   "V8296_021",
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
                                   "D"), check_col_length = F)

# Define stationary periods as breeding, stopover, or non-breeding #############
geo.all <- geo.all %>% group_by(geo_id) %>% mutate(site_type = case_when(
  (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ "Breeding",
  sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "South and partial North" & path.elongated == T ~ "Breeding",
  period == "Non-breeding period" & (duration >= 14 | sitenum == 1 | sitenum == max(sitenum)) ~ "Nonbreeding",
  .default = "Stopover"))

# Calculate IQR between 97.5 and 2.5 quantiles 
geo.all <- geo.all %>% mutate(IQR.lon.dist = distHaversine(cbind(Lon.97.5.,Lat.50.), cbind(Lon.2.5., Lat.50.)),
                              IQR.lat.dist = distHaversine(cbind(Lon.50.,Lat.97.5.), cbind(Lon.50.,Lat.2.5.)))

# Find the time of departure from the breeding grounds #############
geo.all.br.departure <- geo.all %>% filter(period %in% c("Post-breeding migration","Non-breeding period"), 
                                           Recorded_North_South_mig %in% c("Both", "South and partial North", "South")) %>%
  group_by(geo_id) %>% mutate(proximity = distHaversine(cbind(Lon.50., Lat.50.), cbind(deploy.longitude, deploy.latitude))) %>%
  filter(proximity <= 250000) %>% summarize(fall.br.departure = max(EndTime)) 
  #mutate(lon.proximity = deploy.longitude - Lon.50., lat.proximity = deploy.latitude - Lat.50.) %>%
  #filter(abs(lon.proximity) <= 2 & abs(lat.proximity) <= 4) %>% summarize(fall.br.departure = max(EndTime))
  #filter(ifelse(study.site %in% c("Quebec", "Mount Mansfield, Vermont, USA", "Nova Scotia, Canada"), (abs(lon.proximity) < 0.1 & abs(lat.proximity) < 0.1), (abs(lon.proximity) <= 2 & abs(lat.proximity) <= 4))) %>% summarize(fall.br.departure = max(EndTime )) 

# Find the time of return to the breeding grounds #############
geo.all.br.arrival <- geo.all %>%  filter(period %in% c("Pre-breeding migration","Non-breeding period"),
                                          Recorded_North_South_mig %in% c("Both","North"),
                                          path.elongated == F, 
                                          !(geo_id %in% c("V8296_007", "V8296_008"))) %>%
  group_by(geo_id) %>% mutate(proximity = ifelse(geo_id != "WRMA04173", distHaversine(cbind(Lon.50., Lat.50.), cbind(deploy.longitude, deploy.latitude)),
                                                  distHaversine(cbind(Lon.50., Lat.50.), cbind(last(Lon.50.), last(Lat.50.))))) %>%
  filter(proximity <= 250000) %>% summarize(spring.br.arrival = min(StartTime )) 
  #mutate(lon.proximity = ifelse(geo_id != "WRMA04173", deploy.longitude - Lon.50., last(Lon.50.) - Lon.50.),
  #       lat.proximity =ifelse(geo_id != "WRMA04173", deploy.latitude - Lat.50., last(Lat.50.) - Lat.50.)) %>%
  #filter(abs(lon.proximity) <= 2 & abs(lat.proximity) <= 4) %>% summarize(spring.br.arrival = min(StartTime))

# Merge breeding grounds arrival & departure times with location data #######
geo.all <- geo.all %>% dplyr::select(!c(fall.br.departure, spring.br.arrival)) %>% merge(geo.all.br.departure, by = "geo_id", all = T) %>%
  merge(geo.all.br.arrival, by = "geo_id", all = T) %>% group_by(geo_id) %>% 
  arrange(geo_id, StartTime) %>%
  group_by(geo_id) %>%
  mutate(period = ifelse((EndTime <= fall.br.departure | StartTime >= spring.br.arrival) & Recorded_North_South_mig == "Both" & path.elongated == F, "Breeding", period),
         site_type = ifelse((EndTime  <= fall.br.departure | StartTime >= spring.br.arrival) & Recorded_North_South_mig == "Both"& path.elongated == F, "Breeding", site_type),
         period = ifelse(EndTime <= fall.br.departure & (Recorded_North_South_mig %in% c("South and partial North", "South") | path.elongated == T), "Breeding", period),
         site_type = ifelse(EndTime <= fall.br.departure & (Recorded_North_South_mig %in% c("South and partial North", "South") | path.elongated == T), "Breeding", site_type),
         period = ifelse(StartTime >= spring.br.arrival & Recorded_North_South_mig == "North", "Breeding", period),
         site_type = ifelse(StartTime >= spring.br.arrival & Recorded_North_South_mig == "North", "Breeding", site_type))

# For geolocators that were deployed in close proximity, assign one location to the geolocator deployment site
# Here this location will be the median of all geolocator deployment locations for that site. 

# we first add new columns with the median deployment locations in our reference dataset (metadata for each geolcator)
ref_data <- ref_data %>% group_by(study.site) %>% 
  mutate(mod.deploy.lon = median(deploy.longitude))%>%
  mutate(mod.deploy.lat = median(deploy.latitude))

#View(ref_data %>% group_by(study.site, geo.id) %>% summarise(disp = mean(deploy.longitude)))

# Now we can modify our location dataset
#geo.all <- merge(geo.all, ref_data[,c("geo.id", "mod.deploy.lon", "mod.deploy.lat")], by.x = "geo_id", by.y = "geo.id")

geo.all <- geo.all %>% group_by(geo_id) %>% mutate(Lon.50. = case_when(
  (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ mod.deploy.lon,
  sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ mod.deploy.lon,
  #(sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ mod.deploy.lon,
  .default = Lon.50.)) %>%
  mutate(Lat.50. = case_when(
    (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ mod.deploy.lat,
    sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ mod.deploy.lat,
    # (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ mod.deploy.lat,
    .default = Lat.50.))

# Check the location of breeding and nonbreeding sites   
#View(geo.all[,c("geo_id", "Lon.50.", "Lat.50.", "mod.deploy.lon", "mod.deploy.lat")])

# plot locations of deployment sites (there should be 10) 
dep.sites <- geo.all[(geo.all$sitenum == 1),]
plot(wrld_simpl, col = "lightgray", xlim = c(-170, -30), ylim = c(-20, 70))
points(dep.sites$Lon.50., dep.sites$Lat.50., col = "darkred")

# plot locations of all stationary sites 
stat.sites <- geo.all[(geo.all$sitenum > 0),]
plot(wrld_simpl, col = "lightgray", xlim = c(-170, -30), ylim = c(-20, 70))
points(stat.sites$Lon.50., stat.sites$Lat.50., col = "darkred")

# Count stationary sites by season for each geolocator #########################
geo.all <- geo.all %>% arrange(geo_id, StartTime)
num <- 1
geo.all$NB_count <- NA
geo.id.ind <- geo.all$geo_id[1]

for (i in seq(1, nrow(geo.all))){
  
  if (geo.all$site_type[i] == "Nonbreeding" & geo.all$geo_id[i] == geo.id.ind){
    geo.all$NB_count[i] <- num
    num <- num + 1 
  }
  
  if (geo.all$site_type[i] == "Nonbreeding" & geo.all$geo_id[i] != geo.id.ind){
    geo.id.ind <- geo.all$geo_id[i]
    num <-1 
    geo.all$NB_count[i] <- num
    num <- num + 1 
  }
}

# fix dates for geolocators that recorded with the wrong date ##################
geo.all <- geo.all %>% arrange(geo_id, StartTime) %>% 
  group_by(geo_id) %>% mutate(year.difference = year(deploy.on.date) - year(first(StartTime)))

year(geo.all$StartTime) <- year(geo.all$StartTime) + geo.all$year.difference  
year(geo.all$EndTime) <- year(geo.all$EndTime) + geo.all$year.difference  

geo.all$IH.calib.start <- anytime(geo.all$IH.calib.start)
geo.all$IH.calib.end <- anytime(geo.all$IH.calib.end)

year(geo.all$IH.calib.start) <- year(geo.all$IH.calib.start) + geo.all$year.difference  
year(geo.all$IH.calib.end) <- year(geo.all$IH.calib.end) + geo.all$year.difference

geo.all$nbr.arrival <-  as.Date(as.numeric(geo.all$nbr.arrival), optional = T)
geo.all$nbr.departure <- as.Date(as.numeric(geo.all$nbr.departure), optional = T)
geo.all$br.arrival <- as.Date(as.numeric(geo.all$br.arrival), optional = T)
geo.all$br.departure <- as.Date(as.numeric(geo.all$br.departure), optional = T)

year(geo.all$nbr.arrival) <-  year(geo.all$nbr.arrival) + geo.all$year.difference  
year(geo.all$nbr.departure) <- year(geo.all$nbr.departure) + geo.all$year.difference  
year(geo.all$br.arrival) <- year(geo.all$br.arrival) + geo.all$year.difference  
year(geo.all$br.departure) <- year(geo.all$br.departure) + geo.all$year.difference

# Save the data edited here  ###################################################

#Save the new columns in the reference data
write.csv(ref_data,"C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv", row.names=FALSE)

#save the processed geolocatorlocations 
write.csv(geo.all, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/All.locations.csv",row.names=FALSE)

################################################################################
# Fall migration network #######################################################
################################################################################

# Clustering breeding and nonbreeding location #################################

# Either import a file with manual clusters, or create cluster in R

# first we must filter our location data to retain only those relevant to the first breeding period and fall migration 
fall.stat <- geo.all %>% filter(sitenum > 0, duration >= 2, 
                                anydate(StartTime, asUTC = T) <= anydate(nbr.arrival),
                                Recorded_North_South_mig %in% c("Both", "South and partial North", "South"))

# we will also extract all fall locations (including while the bird was moving), for figures
fall.move <- geo.all %>% filter(site_type %in% anydate(StartTime, asUTC = T) <= anydate(nbr.arrival),
                                Recorded_North_South_mig %in% c("Both", "South and partial North", "South"))

#get the timing of the first nonbreeding area
fall.timings.nb <- geo.all %>% filter(NB_count == 1) %>% dplyr::select(NB.first.site.arrival = StartTime)
fall.stat <- merge(fall.stat, fall.timings.nb, by = "geo_id")

#only retain the stopovers and the first nonbreeding sites occupied, and filter out stopovers that are within 250 km of breeding site 
fall.stat <- fall.stat %>% group_by(geo_id) %>% filter(StartTime <= NB.first.site.arrival) %>% 
  # filter(case_when(!(study.site %in% c("Quebec", "Mount Mansfield, Vermont, USA", "Nova Scotia, Canada")) ~ site_type != "Breeding",
  #                  T ~ site_type == site_type )) %>%
  # filter(ifelse(study.site %in% c("Quebec", "Mount Mansfield, Vermont, USA", "Nova Scotia, Canada"), distHaversine(cbind(Lon.50.,Lat.50.), cbind(deploy.longitude, deploy.latitude)) > 50000, distHaversine(cbind(Lon.50.,Lat.50.), cbind(deploy.longitude, deploy.latitude)) > 250000)) %>%
  # mutate(site_type = ifelse(site_type == "Breeding", "Stopover", site_type),
  #        period  = ifelse(site_type == "Breeding", "Post-breeding migration", period),)
  filter(distHaversine(cbind(Lon.50.,Lat.50.), cbind(deploy.longitude, deploy.latitude)) > 250000) %>%
  mutate(site_type = ifelse(site_type == "Breeding", "Stopover", site_type),
         period = ifelse(period == "Breeding", "Post-breeding migration", period))

  
# Create clusters in two steps to account for the equinox 

# load equinox polygon 
equipol <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Manual_stat_site_clustering/Layers/equipol.shp", crs = st_crs(wrld_simpl))
fall.stat.sf <- st_as_sf(fall.stat, coords = c("Lon.50.", "Lat.50."), crs = st_crs(wrld_simpl), remove = F)

# separate points affected by the equinox and those that are not 
fall.stat.equi <- st_intersection(fall.stat.sf, equipol)
fall.stat.equi <- st_drop_geometry(fall.stat.equi)

fall.stat.norm <- st_difference(fall.stat.sf, equipol)
fall.stat.norm <- st_drop_geometry(fall.stat.norm )

# cluster points in each group separately, then merge the cluster info 
cluster.data1 <- clusterLocs(locs = fall.stat.equi, maxdiam = 743, lon.only = T)
cluster.data2 <- clusterLocs(locs = fall.stat.norm, maxdiam = 700)

cluster.data2$clusters <- cluster.data2$clusters + max(cluster.data1$clusters)

fall.stat.equi$cluster <- cluster.data1$clusters
fall.stat.norm$cluster <- cluster.data2$clusters

fall.stat <- rbind(fall.stat.equi, fall.stat.norm)

# # Uncomment this code to generate clusters using the pam function
# cluster.data <- clusterLocs(locs = fall.stat, maxdiam = 700)
# fall.stat$cluster <- cluster.data$clusters
# 
# # Export fall stopovers for manual clustering in QGIS
# fall.stat.sites <- st_as_sf(data.frame(fall.stat[,c(1:12, which(colnames(fall.stat) == "cluster"))]), coords = c("Lon.50.", "Lat.50."))
# st_crs(fall.stat.sites) <- st_crs(wrld_simpl)
# st_write(fall.stat.sites, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Manual_stat_site_clustering/layers/fall_stat_sites6.shp", append=FALSE)

# # Import clusters created manually
# fall.manual.cluster <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Manual_stat_site_clustering/Tables/Fall_manual_clusters_conservativeV6.csv") %>%
#   rename(cluster = Cluster, cluster.region = ClusterReg) %>%
#   mutate(cluster.region= as.factor(cluster.region)) %>%
#   mutate(cluster = as.numeric(cluster.region))

# # Merge manual cluster info with original dataset
# fall.stat <- merge(fall.stat, fall.manual.cluster[,c("geo_id", "sitenum", "cluster", "cluster.region")], by = c("geo_id", "sitenum"))

# Add breeding sites with separate clusters
fall.breed.n <- geo.all  %>% group_by(geo_id) %>% filter(StartTime == first(StartTime), geo_id %in% unique(fall.stat$geo_id)) %>% 
  group_by(study.site) %>% summarize(geo_per_site = n())

fall.breed <- geo.all  %>% filter(StartTime == first(StartTime), geo_id %in% unique(fall.stat$geo_id)) %>%
  arrange(study.site, geo_id, StartTime)

fall.breed$cluster <- rep(seq(max(fall.stat$cluster) + 1, max(fall.stat$cluster) + nrow(fall.breed.n)), fall.breed.n$geo_per_site)

fall.stat <- bind_rows(fall.stat, fall.breed) %>% arrange(geo_id, StartTime)

# plot stopover and nonbreeding nodes
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = site_type)) +
  #geom_text(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., label = geo_id), cex = 2)+
  theme_bw() +
  theme(text = element_text(size = 16))

# plot the fall clusters 
reg.bounds <- st_read("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/Spatial_data_downloads/ne_50m_admin_1_states_provinces/ne_50m_admin_1_states_provinces.shp")

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  geom_sf(data = reg.bounds, fill = NA, lwd = 0.2, alpha = 1) +
   coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = fall.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  geom_errorbar(data = fall.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  #geom_path(data = fall.stat[fall.stat$geo_id == "A",], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  #geom_point(data = fall.stat[fall.stat$geo_id == "A",], mapping = aes(x = Lon.50., y = Lat.50.), alpha = 0.5) +
  #geom_text(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., label = geo_id))+
  geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5, linewidth = 0.1) +
  geom_point(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 2) +
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None")

# Optional: add edge from nonbreeding site back to the breeding site of origin
fall.breed.return <- fall.breed %>% group_by(geo_id) %>%
  mutate(sitenum = max(fall.stat[(fall.stat$geo_id == geo_id),]$sitenum))
fall.breed.return[,c("StartTime", "EndTime", "duration")] <- NA
fall.breed.return[,c("period")] <- "Post-breeding migration"
fall.stat <- bind_rows(fall.stat, fall.breed.return) %>% arrange(geo_id, sitenum)

# Generate the network from our location data and clusters #####################

# Add a column with inter-cluster connections to our location dataset
fall.stat <- fall.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# Create a dataframe with the edge info 
fall.edge.df <- fall.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                            Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                            Lat.2.5., Lat.97.5., EndTime,
                                            sitenum, duration, period,study.site,
                                            Range_region, NB_count, period,
                                            site_type) %>%
  filter(!is.na(next.cluster))

# calculate node locations as the median locations of the points included within 
# Note, we use the 50% quantile of the posterior distribution 
fall.node.lon <- fall.stat %>% group_by(cluster) %>% 
  summarize(node.lon = median(Lon.50.))

fall.node.lat <- fall.stat %>% group_by(cluster) %>% 
  summarize(node.lat = median(Lat.50.))

# We also assess which use was the more common for a node (breeding, nonbreeding, or stopover)
# The most common use will be the one assigned to the node
# fall.node.type <- fall.stat %>% group_by(cluster, site_type) %>%
#   summarize(use = n()) %>% filter(use == max(use)) %>% 
#   mutate(site_type = factor(site_type, levels = c("Breeding", "Nonbreeding", "Stopover")))%>%
#   arrange(site_type) %>%
#   mutate(site_type_num = case_when(
#     site_type == "Breeding" ~ 1,
#     site_type == "Nonbreeding" ~ 2,
#     site_type == "Stopover" ~ 3
#   )) %>% distinct(cluster, .keep_all = TRUE) %>% arrange(cluster) #must verify whether sites have equal maxima of uses

fall.node.type <- fall.stat %>% group_by(cluster) %>%
  summarize(count.breeding = length(site_type[site_type == "Breeding"]),
            count.stopover = length(site_type[site_type == "Stopover"]),
            count.nonbreeding = length(site_type[site_type == "Nonbreeding"])) %>%
  mutate(site_type_num = case_when(
    count.breeding > 0  ~ 1,
    count.stopover > 0 & count.nonbreeding == 0 ~ 2,
    count.nonbreeding > 0 ~ 3)) %>%
  mutate(site_type = case_when(
    site_type_num == 1  ~ "Breeding",
    site_type_num == 2 ~ "Stopover",
    site_type_num == 3 ~ "Nonbreeding")) %>%
  mutate(site_type = factor(site_type, levels = c("Breeding", "Nonbreeding", "Stopover")))

# Create a layout with node locations 
meta.fall <- data.frame("vertex" = seq(1, max(fall.edge.df$cluster)), 
                        #data.frame("vertex" = seq(min(fall.edge.df$cluster), max(fall.edge.df$cluster)), 
                        "Lon.50." = fall.node.lon$node.lon,
                        "Lat.50." = fall.node.lat$node.lat,
                        "node.type" = fall.node.type$site_type,
                        "node.type.num" = fall.node.type$site_type_num)

# For fall nodes where latitudinal accuracy is low, set location close to the coast
meta.fall[c(3, 5, 7),]$Lat.50. <- c(35.3, 41.45, 44.77)

fall.location <- as.matrix(meta.fall[, c("Lon.50.", "Lat.50.")])

# Create the fall network
fall.graph <- graph_from_data_frame(fall.edge.df, directed = T, vertices = meta.fall)

# Colour palette for site type
type.palette <- rainbow(3)  

# plot the fall network over North and South America
plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph, vertex.size = 200, vertex.size2 = 200,
     edge.width = 1, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.color = type.palette[meta.fall$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.fall$node.type.num)],
       pch = 16)

# Add edge weights to the fall network to show the flow of individuals #########

# We start with edges based on the number of individuals moving between nodes

# list of connections and the number of times they occur
fall.con <- fall.edge.df %>% group_by(cluster, next.cluster) %>% 
  summarize(weight = n()) 

# Create a fall network with weighed edges
fall.graph.weighed <- graph_from_data_frame(fall.con, directed = T, vertices = meta.fall)
is_weighted(fall.graph.weighed)

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con$weight/1.5, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con)),
     vertex.color = type.palette[meta.fall$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.fall$node.type.num)],
       pch = 16)

################################################################################
# Spring migration network #####################################################
################################################################################

# Clustering breeding and nonbreeding location #################################

# Either import a file with manual clusters, or create cluster in R

# Create clusters for spring stopover locations and last nonbreeding locations in R 
spring.stat <- geo.all %>% filter(sitenum > 0, duration > 2,# site_type %in% c("Stopover","Nonbreeding"),
                                  period %in% c("Pre-breeding migration","Non-breeding period"),
                                  Recorded_North_South_mig %in% c("Both","North", "South and partial North"),
                                  !(geo_id %in% c("V8296_007", "V8296_008")))

# we will also extract all fall mlocations (including while the bird was moving), for figures
spring.move <- geo.all %>% filter(site_type %in% c("Stopover","Nonbreeding"),
                                  period %in% c("Pre-breeding migration","Non-breeding period"),
                                  Recorded_North_South_mig %in% c("Both","North", "South and partial North"),
                                  !(geo_id %in% c("V8296_007", "V8296_008")))

#get the timing of the last nonbreeding area
spring.timings.nb <- geo.all %>% group_by(geo_id) %>% filter(NB_count == max(NB_count, na.rm = T)) %>% dplyr::select(NB.last.site.arrival = StartTime)
spring.stat <- merge(spring.stat, spring.timings.nb, by = "geo_id")

#only retain the stopovers and the first nonbreeding sites occupied, and filter out stopovers that are within 250 km of breeding site 
spring.stat <- spring.stat %>% group_by(geo_id) %>% filter(StartTime >= NB.last.site.arrival) %>% group_by(geo_id)%>% 
  filter(distHaversine(cbind(Lon.50.,Lat.50.), cbind(deploy.longitude, deploy.latitude)) > 250000 | geo_id == "WRMA04173")%>%
  mutate(site_type = ifelse(site_type == "Breeding", "Stopover", site_type),
         period = ifelse(period == "Breeding", "Post-breeding migration", period))

# Uncomment this code to generate clusters using the pam function
cluster.data <- clusterLocs(locs = spring.stat, maxdiam = 732)
spring.stat$cluster <- cluster.data$clusters

# # export spring stat sites for manual clustering
# spring.stat.sites <- st_as_sf(spring.stat[,c(1:12, which(colnames(spring.stat) == "cluster"))], coords = c("Lon.50.", "Lat.50."))
# st_crs(spring.stat.sites) <- st_crs(wrld_simpl)
# st_write(spring.stat.sites, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Manual_stat_site_clustering/Layers/spring_stat_sitesV3.shp", append = F)

# # Import clusters created manually
# spring.manual.cluster <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Manual_stat_site_clustering/Tables/Spring_manual_clusters_V2.csv")
# spring.manual.cluster <- spring.manual.cluster %>% rename(cluster = Cluster, cluster.region = ClusterReg) %>%
#   mutate(cluster.region= as.factor(cluster.region)) %>%
#   mutate(cluster = as.numeric(cluster.region))
# 
# # Merge manual cluster info with original dataset
# spring.stat <- merge(spring.stat, spring.manual.cluster[,c("geo_id", "sitenum", "cluster", "cluster.region")], by = c("geo_id", "sitenum"))

# Add breeding sites with separate clusters
spring.breed.n <- geo.all  %>% group_by(geo_id) %>% filter(geo_id %in% unique(spring.stat$geo_id), StartTime == first(StartTime)) %>% 
  group_by(study.site) %>% summarize(geo_per_site = n())

spring.breed <- geo.all %>% filter(geo_id %in% unique(spring.stat$geo_id), site_type == "Breeding",
                                   StartTime == last(StartTime)) %>%
  arrange(study.site, geo_id, StartTime)

spring.breed$cluster <- rep(seq(max(spring.stat$cluster) + 1, max(spring.stat$cluster) + nrow(spring.breed.n)), spring.breed.n$geo_per_site)

spring.stat <- bind_rows(spring.stat, spring.breed) %>% arrange(geo_id, StartTime)

# plot stationary stopover and nonbreeding sites 
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  geom_sf(data = reg.bounds, fill = NA, lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = Breeding_region_MC))+ 
  theme_bw() +
  theme(text = element_text(size = 14))

# plot the clusters 
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  geom_sf(data = reg.bounds, fill = NA, lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = spring.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  geom_errorbar(data = spring.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  #geom_path(data = spring.stat[spring.stat$geo_id == "4105_008",], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  #geom_point(data = spring.stat[spring.stat$geo_id == "4105_008",], mapping = aes(x = Lon.50., y = Lat.50.), alpha = 0.5) +
  #geom_text(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., label = geo_id), cex = 4, nudge_x =  1)+
  geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5, linewidth = 0.1) +
  geom_point(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 1.5) +
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None")

# Optional: add edge from nonbreeding site back to the breeding site of origin
spring.breed.return <- spring.breed %>% group_by(geo_id) %>%
  mutate(sitenum = 1)
spring.breed.return[,c("StartTime", "EndTime", "duration")] <- NA
spring.breed.return[,c("period")] <- "Pre-breeding migration"

# For geolocators deployed in the nonbreeding grounds, we need to increase site numbers because the breeding site
spring.stat <- spring.stat %>% group_by(geo_id) %>% mutate(sitenum = ifelse(deploy.range == "Nonbreeding", sitenum + 1, sitenum))
spring.stat <- bind_rows(spring.stat, spring.breed.return) %>% arrange(geo_id, sitenum)

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = spring.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  geom_errorbar(data = spring.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, alpha = 0.3, color = "black") +
  geom_point(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 1.5) +
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None")

# Generate the network from our location data and clusters #####################

# Add a column with inter-cluster connections to our location dataset
spring.stat  <- spring.stat %>% mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# Create a dataframe with the edge info 
spring.edge.df <- spring.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                                Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                                Lat.2.5., Lat.97.5., EndTime,
                                                sitenum, duration, period,study.site,
                                                Range_region, NB_count, period,
                                                site_type) %>%
  filter(!is.na(next.cluster))

# calculate node locations as the median locations of the points included within 
# Note, we use the 50% quantile of the posterior distribution 
spring.node.lon <- spring.stat %>% group_by(cluster) %>% 
  summarize(node.lon = median(Lon.50.))

spring.node.lat <- spring.stat %>% group_by(cluster) %>% 
  summarize(node.lat = median(Lat.50.))

# We also assess which use was the more common for a node (breeding, nonbreeding, or stopover)
# The most common use will be the one assigned to the node
# spring.node.type <- spring.stat %>% group_by(cluster, site_type) %>%
#   summarize(use = n()) %>% filter(use == max(use)) %>%
#   mutate(site_type = factor(site_type, levels = c("Breeding", "Nonbreeding", "Stopover")))%>%
#   arrange(site_type) %>%
#   mutate(site_type_num = case_when(
#     site_type == "Breeding" ~ 1,
#     site_type == "Nonbreeding" ~ 2,
#     site_type == "Stopover" ~ 3
#   )) %>% distinct(cluster, .keep_all = TRUE) %>% arrange(cluster) #must verify whether sites have equal maxima of uses

spring.node.type <- spring.stat %>% group_by(cluster) %>%
  summarize(count.breeding = length(site_type[site_type == "Breeding"]),
            count.stopover = length(site_type[site_type == "Stopover"]),
            count.nonbreeding = length(site_type[site_type == "Nonbreeding"])) %>%
  mutate(site_type_num = case_when(
    count.breeding > 0  ~ 1,
    count.stopover > 0 & count.nonbreeding == 0 ~ 2,
    count.nonbreeding > 0 ~ 3)) %>%
  mutate(site_type = case_when(
    site_type_num == 1  ~ "Breeding",
    site_type_num == 2 ~ "Stopover",
    site_type_num == 3 ~ "Nonbreeding")) %>%
  mutate(site_type = factor(site_type, levels = c("Breeding", "Nonbreeding", "Stopover")))

# Create a layout with node locations 
meta.spring <- data.frame("vertex" = seq(min(spring.stat$cluster), max(spring.stat$cluster)), 
                          "Lon.50." = spring.node.lon$node.lon,
                          "Lat.50." = spring.node.lat$node.lat,
                          "node.type.num" = spring.node.type$site_type_num,
                          "node.type" = spring.node.type$site_type)

spring.location <- as.matrix(meta.spring[, c("Lon.50.", "Lat.50.")])

# Create the spring network
spring.graph <- graph_from_data_frame(spring.edge.df, directed = T, vertices = meta.spring)

# Colour palette for site type
type.palette <- rainbow(3)  

# plot the spring network over North and South America
plot(wrld_simpl, xlim = c(-170, -30),ylim = c(-15, 70))
plot(spring.graph, vertex.size = 200, vertex.size2 = 200,
     edge.width = 1, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.color = type.palette[meta.spring$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.spring$node.type.num)],
       pch = 16)

# Add edge weights to the spring network to show the flow of individuals #########

# We start with edges based on the number of individuals moving between nodes

# list of connections and the number of times they occur
spring.con <- spring.edge.df %>% group_by(cluster, next.cluster) %>% 
  summarize(weight = n()) 

spring.con[spring.con$cluster == 25 | spring.con$next.cluster == 25,]

# Create a new spring network
spring.graph.weighed <- graph_from_data_frame(spring.con, directed = T, vertices = meta.spring)
is_weighted(spring.graph.weighed)

plot(wrld_simpl, xlim = c(-170, -30),ylim = c(-15, 70))
plot(spring.graph.weighed, vertex.size = 200, vertex.size2 = 200,vertex.label.dist = 30,
     edge.width = spring.con$weight/1.5, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.3, 0.3), nrow(spring.con)),
     vertex.color = type.palette[meta.spring$node.type.num], add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.spring$node.type.num)],
       pch = 16)

# ################################################################################
# # Nonbreeding movements network ################################################
# ################################################################################
# 
# # Clustering breeding and nonbreeding location #################################
# 
# # Either import a file with manual clusters, or create cluster in R
# 
# # Create clusters for nonbreeding sites in R 
# NB.stat <- geo.all %>% filter(sitenum > 0, site_type %in% c("Nonbreeding"),
#                                   period %in% c("Non-breeding period"),
#                                   Recorded_North_South_mig %in% c("Both", "South and partial North"))
# 
# # #get the timing of the first nonbreeding area
# # NB.start.time <- geo.all %>% filter(NB_count == 1) %>% dplyr::select(NB.first.site.arrival = StartTime)
# # NB.stat <- merge(NB.stat, NB.start.time, by = "geo_id")
# # 
# # #get the timing of the lastnonbreeding area
# # NB.end.time <- geo.all %>% group_by(geo_id) %>% filter(NB_count == max(NB_count, na.rm = T)) %>% dplyr::select(NB.last.site.arrival = StartTime)
# # NB.stat <- merge(NB.stat.stat, NB.end.time.nb, by = "geo_id")
# 
# # code by stack overflow user jlhoward (https://stackoverflow.com/questions/21095643/approaches-for-spatial-geodesic-latitude-longitude-clustering-in-r-with-geodesic)
# # Calculates the geodesic distance between points and creates a distance matrix
# geo.dist = function(df) {
#   require(geosphere)
#   d <- function(i,z){         # z[1:2] contain long, lat
#     dist <- rep(0,nrow(z))
#     dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
#     return(dist)
#   }
#   dm <- do.call(cbind,lapply(1:nrow(df),d,df))
#   return(as.dist(dm))
# }
# 
# dist.matrix <- geo.dist(NB.stat[,c("Lon.50.","Lat.50.")])
# NB.clust <- hclust(dist.matrix)
# plot(NB.clust)
# NB.stat$cluster  <- cutree(NB.clust, k = 5)
# 
# # plot stopover and nonbreeding nodes 
# 
# ggplot(st_as_sf(wrld_simpl))+
#   geom_sf(colour = NA, fill = "lightgray") +
#   coord_sf(xlim = c(-90, -30),ylim = c(-15, 20))+
#   geom_errorbar(data = NB.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), color = "red", width=1, alpha = 0.2) + 
#   geom_point(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = site_type)) +
#   geom_path(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5,
#             arrow = arrow(end = "last", type = "open", length = unit(0.10, "inches"))) 
# 
# ggplot(st_as_sf(wrld_simpl))+
#   geom_sf(colour = NA, fill = "lightgray") +
#   coord_sf(xlim = c(-90, -30),ylim = c(-15, 20)) +
#   geom_point(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = cluster)) +
#   geom_path(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5,
#             arrow = arrow(end = "last", type = "open", length = unit(0.10, "inches"))) +
#    scale_colour_gradientn(colours=rainbow(5))
# 
# 
# ggplot(st_as_sf(wrld_simpl))+
#   geom_sf(colour = NA, fill = "lightgray") +
#   coord_sf(xlim = c(-90, -30),ylim = c(-15, 20)) +
#   geom_point(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = Range_region)) +
#   geom_path(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5,
#             arrow = arrow(end = "last", type = "open", length = unit(0.10, "inches"))) 
# 
# # Generate the network from our location data and clusters #####################
# 
# # Add a column with inter-cluster connections to our location dataset
# NB.stat  <- NB.stat %>% mutate(next.cluster = case_when(
#   lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
#   .default = NA)) 
# 
# # Create a dataframe with the edge info 
# NB.edge.df <- NB.stat %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
#                                                 Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
#                                                 Lat.2.5., Lat.97.5., EndTime,
#                                                 sitenum, duration, period,study.site,
#                                                 Range_region, NB_count, period,
#                                                 site_type) %>%
#   filter(!is.na(next.cluster))
# 
# # calculate node locations as the median locations of the points included within 
# # Note, we use the 50% quantile of the posterior distribution 
# NB.node.lon <- NB.stat %>% group_by(cluster) %>% 
#   summarize(node.lon = median(Lon.50.))
# 
# NB.node.lat <- NB.stat %>% group_by(cluster) %>% 
#   summarize(node.lat = median(Lat.50.))
# 
# # We also assess which use was the more common for a node (breeding, nonbreeding, or stopover)
# # The most common use will be the one assigned to the node
# NB.node.type <- NB.stat %>% group_by(cluster, site_type) %>%
#   summarize(use = n()) %>% filter(use == max(use)) %>% mutate(site_type_num = case_when(
#     site_type == "Nonbreeding" ~ 1,
#     site_type == "Stopover" ~ 2
#   ))
# 
# # Create a layout with node locations 
# meta.NB <- data.frame("vertex" = seq(min(NB.stat$cluster), max(NB.stat$cluster)), 
#                    "Lon.50." = NB.node.lon$node.lon,
#                    "Lat.50." = NB.node.lat$node.lat,
#                    "node.type.num" = NB.node.type$site_type_num,
#                    "node.type" = NB.node.type$site_type)
# 
# NB.location <- as.matrix(meta.NB[, c("Lon.50.", "Lat.50.")])
# 
# # Create the spring network
# NB.graph <- graph_from_data_frame(NB.edge.df, directed = T, vertices = meta.NB)
# 
# # Colour palette for site type
# type.palette <- rainbow(3)  
# 
# # plot the fall network over North and South America
# plot(NB.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
#      edge.width = 1, edge.arrow.size = 0.5, edge.arrow.width = 0.5,  
#      layout = NB.location, rescale = F, asp = 0, xlim = c(-90, -30),
#      ylim = c(-15, 20), vertex.color = type.palette[meta.NB $node.type.num])
# plot(wrld_simpl, add = T)
# 
# # Add edge weights to the spring network to show the flow of individuals #########
# 
# # We start with edges based on the number of individuals moving between nodes
# 
# # list of connections and the number of times they occur
# NB.con <- NB.edge.df %>% group_by(cluster, next.cluster) %>% 
#   summarize(weight = n()) 
# 
# # Create a  network for nonbreeding movements 
# NB.graph.weighed <- graph_from_data_frame(NB.con, directed = T, vertices = meta.NB )
# is_weighted(NB.graph.weighed)
# 
# plot(wrld_simpl, xlim = c(-90, -30), ylim = c(-15, 20))
# plot(NB.graph.weighed, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
#      edge.width = NB.con$weight, edge.arrow.size = 0.10, edge.arrow.width = 1,  
#      layout = NB.location, rescale = F, asp = 0, edge.curved = rep(c(-0.3, 0.3), nrow(NB.con)),
#      vertex.color = type.palette[meta.NB $node.type.num], add = T)

################################################################################
# Assess and propagate the relative abundance of blackpoll warblers from eBird
################################################################################

# create and export geolocator deployment sites
geo.breed <- geo.all %>% filter(site_type == "Breeding", (StartTime == first(StartTime) | StartTime == last(StartTime))) %>%
  distinct(geo_id, .keep_all = T) %>%
  dplyr::select(geo_id, Lon.50., Lat.50., sitenum, site_type, period, study.site)
breed.sites <- st_as_sf(geo.breed, coords = c("Lon.50.", "Lat.50."))

# set crs of breeding.sites sf object
st_crs(breed.sites) <- st_crs(wrld_simpl)

# export breeding sites for to draw breeding region polygons.  
# st_write(breed.sites, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_breeding_sites.shp")

# import breeding regions polygon and join attributes with breeding site data 
abundance.regions <- read_sf("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/Relative_abundance_propagation/bpw_abundance_regions_adjusted.shp")
abundance.regions <- st_join(abundance.regions, breed.sites)

# Add the abundance region of each individual to the reference dataset. These will be the same regions used to calculate migratory connectivity 
ref_data <- ref_data %>% dplyr::select(-Breeding_region_MC)
ref_data <- merge(ref_data, data.frame(abundance.regions)[,c("region", "geo_id")], by.x = "geo.id", by.y = "geo_id" ) %>%
  rename(Breeding_region_MC = region)

write.csv(ref_data,"C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv", row.names=FALSE)

# import breeding season abundance data
bpw.fall.ab <- load_raster(path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/geo_spatial_data/eBird_imports",
                           species = "bkpwar",
                           product = "abundance",
                           period = "seasonal",
                           resolution = "9km")

bpw.fall.ab <- terra::project(bpw.fall.ab, crs(abundance.regions))

# plot breeding regions, breeding sites/geolocator deployment sites, and abundance 
# NOTE: the breeding site for WRMA04173 in Northern Quebec will only be used for the spring migration
abundance.regions <- st_transform(abundance.regions , crs = crs(wrld_simpl))
plot(bpw.fall.ab$breeding, xlim = c(-180, -30), ylim = c(30, 80), legend = FALSE)
plot(as_Spatial(abundance.regions), col = NA, border = "darkred", lwd = 3, add = T)
plot(wrld_simpl[wrld_simpl$NAME %in% c("United States", "Canada"),], add = T)
#points(geo.breed[geo.breed$geo_id != "WRMA04173",]$Lon.50., geo.breed[geo.breed$geo_id != "WRMA04173",]$Lat.50., cex = 1, col = "black", pch = 19)
points(geo.breed$Lon.50., geo.breed$Lat.50., cex = 1, col = "black", pch = 19)

# extract the abundance for each region,
ab.extract <- terra::extract(bpw.fall.ab$breeding, abundance.regions, fun = sum, na.rm=TRUE)
ab.extract$ID <- abundance.regions$geo_id
ab.extract$breedregionname <- abundance.regions$region

# Create a dataframe with the relative abundance per region
ab.per.region <- merge(as.data.frame(abundance.regions), ab.extract, by.x = "geo_id", by.y = "ID") %>%
  dplyr::select(-geometry, -region) %>%
  rename(br.region.r.abundance = breeding, br.polygon = breedregionname) %>%
  mutate(br.region.prop.total.population = br.region.r.abundance/ sum(unique(br.region.r.abundance)))

################################################################################
# relative abundance propagation on in the fall season 
################################################################################

# calculate the relative abundance for each breeding node (relative abundance of region/number of breeding nodes in region)
# then create abundance units for each individual bird tracked in the fall
# these units will be used to propagate the relative abundance from the breeding regions
fall.breed.ab <- merge(fall.breed, dplyr::select(ab.per.region, geo_id, br.region.prop.total.population, br.polygon), by = "geo_id") %>%
  group_by(br.polygon) %>% mutate(ab.unit = br.region.prop.total.population/n_distinct(geo_id))

#View(fall.breed.ab %>% dplyr::select(geo_id, StartTime, EndTime, Lon.50.,
#                                      Lat.50., br.polygon,
#                                      br.region.prop.total.population,
#                                      br.node.prop.total.population,
#                                      ab.unit))

# We add the abundance units to our dataset of movements
fall.edge.df.ab <- merge(fall.edge.df, dplyr::select(fall.breed.ab, geo_id, ab.unit, br.region.prop.total.population, br.polygon))

# list of connections weighed by abundance unit
# We add an edge type if optional spring edges were added to distinguish fall and spring edges 
fall.con.ab <- fall.edge.df.ab %>% group_by(geo_id) %>%
  mutate(edge.type = if_else(next.cluster == first(cluster), "spring", "fall")) %>%
  group_by(cluster, next.cluster) %>% 
  reframe(weight = sum(ab.unit), weight.n = n(),  edge.type = first(edge.type))

# weight: weight based on abundance units 
# weight.n: weight based on number of tracked individuals moving between nodes 

# Fall graph weighed using eBird relative abundance 

# Create a fall network with weighed edges
fall.graph.weighed.ab <- graph_from_data_frame(fall.con.ab, directed = T, vertices = meta.fall)

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(fall.graph.weighed.ab, vertex.size = 200, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall$node.type.num], vertex.label.dist = 30,
     add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.fall$node.type.num)],
       pch = 16)

################################################################################
# relative abundance propagation on in the spring season 
################################################################################

# calculate the relative abundance for each breeding node (relative abundance of region/number of breeding nodes in region)
# then create abundance units for each individual bird tracked in the spring
# these units will be used to propagate the relative abundance from the breeding regions
spring.breed.ab <- merge(spring.breed, dplyr::select(ab.per.region, geo_id, br.region.prop.total.population, br.polygon), by = "geo_id") %>%
  group_by(br.polygon) %>% mutate(ab.unit = br.region.prop.total.population/n_distinct(geo_id))

# View(spring.breed.ab %>% dplyr::select(geo_id, StartTime, EndTime, Lon.50.,
#                                      Lat.50., br.polygon,
#                                      br.region.prop.total.population,
#                                      br.node.prop.total.population,
#                                      ab.unit))

# We add the abundance units to our dataset of movements
spring.edge.df.ab <- merge(spring.edge.df, dplyr::select(spring.breed.ab, geo_id, ab.unit, br.region.prop.total.population, br.polygon))

# list of connections weighed by abundance unit
# We add an edge type if optional fall edges were included to distinguish fall and spring edges 
spring.con.ab <- spring.edge.df.ab %>% group_by(geo_id) %>%
  mutate(edge.type = if_else(cluster == last(next.cluster), "fall", "spring")) %>%
  group_by(cluster, next.cluster) %>% 
  reframe(weight = sum(ab.unit), weight.n = n(), edge.type = first(edge.type))

# weight: weight based on abundance units 
# weight.n: weight based on number of tracked individuals moving between nodes 

# spring graph weighed using eBird relative abundance 

# Create a spring network with weighed edges
spring.graph.weighed.ab <- graph_from_data_frame(spring.con.ab, directed = T, vertices = meta.spring)

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph.weighed.ab, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*40, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = type.palette[meta.spring$node.type.num], vertex.label.dist = 30, 
     add = T)
legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
       col = type.palette[unique(meta.spring$node.type.num)],
       pch = 16)

################################################################################
# Population proportion by site during the fall season 
################################################################################

# Create a dataframe with the abundance units of individuals from each section of the blackpoll warbler's range in each node
fall.stat.ab <- fall.stat %>% filter(!(is.na(next.cluster))) %>% merge(fall.breed.ab[,c("geo_id", "ab.unit")], by = "geo_id") 

# calculate proportion of abundance units from each region at each node  
fall.ab.by.origin <- fall.stat.ab %>% group_by(cluster, Breeding_region_MC) %>%
  summarize(region.ab.units = sum(ab.unit), region.n = n_distinct(geo_id)) %>% ungroup() %>%
  complete(cluster, Breeding_region_MC, fill = list(region.ab.units = 0)) %>%
  complete(cluster, Breeding_region_MC, fill = list(region.n = 0))

# Convert data from wide to long
fall.ab.by.origin <- fall.ab.by.origin %>% pivot_wider(names_from = Breeding_region_MC, values_from = c(region.ab.units, region.n)) %>%
  rename(prop.ab.central = `region.ab.units_Central Region`,
         prop.ab.eastern = `region.ab.units_Eastern Region`,
         prop.ab.western = `region.ab.units_Western Region`,
         prop.ab.northwestern = `region.ab.units_Northwestern Region`,
         n.central = `region.n_Central Region`,
         n.eastern = `region.n_Eastern Region`,
         n.western = `region.n_Western Region`,
         n.northwestern = `region.n_Northwestern Region`)

# Merge abundance data with the node metadata
meta.fall.ab <- merge(meta.fall, fall.ab.by.origin, by.x = "vertex", by.y = "cluster") 

# Also add a column with the overall abundance at each site, and the number of individual at each site 
fall.stat.ab.per.site <- fall.stat.ab %>% group_by(cluster) %>% 
  summarize(r.abundance.at.cluster = sum(ab.unit))

meta.fall.ab <- merge(meta.fall.ab, fall.stat.ab.per.site, by.x = "vertex", by.y = "cluster") %>% 
  rowwise() %>%
  mutate(n.individuals.at.cluster = sum(n.central, n.eastern, n.western, n.northwestern))  

#create a column that can be converted to a numeric vector
meta.fall.ab <- transform(meta.fall.ab, num.reg.ab.vector = asplit(cbind(prop.ab.central, prop.ab.eastern, prop.ab.western, prop.ab.northwestern), 1))

# Plot of proportional node use during the fall migration 

# Create vector of vertex shapes 
meta.fall.ab <- meta.fall.ab %>% rowwise() %>%
  mutate(shape_single = length(which(c(prop.ab.central, 
                                       prop.ab.eastern, 
                                       prop.ab.western, 
                                       prop.ab.northwestern)==0))) %>%
  ungroup() %>%
  mutate(shape_single = ifelse(shape_single == 3 & node.type != "Breeding", "circle", "none")) %>%
  mutate(shape_single_breeding = ifelse(node.type == "Breeding", "square", "none"))%>%
  mutate(shape_multiple = ifelse(shape_single == "none" & shape_single_breeding == "none", "pie", "none")) %>%
  mutate(shape_colour_single = case_when(shape_single != "none" & prop.ab.central != 0 ~ "#0072B2",
                                         shape_single != "none" & prop.ab.eastern != 0 ~ "#009E73",
                                         shape_single != "none" & prop.ab.western != 0 ~ "#D55E00",
                                         shape_single != "none" & prop.ab.northwestern != 0 ~"#F0E442",
                                         .default = NA)) %>%
  mutate(shape_colour_single_breeding = case_when(shape_single_breeding != "none" & prop.ab.central != 0 ~ "#0072B2",
                                                  shape_single_breeding != "none" & prop.ab.eastern != 0 ~ "#009E73",
                                                  shape_single_breeding != "none" & prop.ab.western != 0 ~ "#D55E00",
                                                  shape_single_breeding != "none" & prop.ab.northwestern != 0 ~ "#F0E442",
                                                  .default = NA))

# Create a palette for site use by range region  
reg.ab.palette <- list(c("#0072B2", "#009E73", "#D55E00", "#F0E442"))

# Prepare edge colours for spring and fall edges (spring edges should not appear in this plot)
edge.cols.fall <- fall.con.ab %>% mutate(col = case_when(
  edge.type == "fall" ~ adjustcolor("darkgray", alpha.f = 0.9),
  edge.type == "spring" ~ NA))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

plot(fall.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
     vertex.shape = meta.fall.ab$shape_single_breeding, vertex.color = meta.fall.ab$shape_colour_single_breeding,
     edge.arrow.width = 0,edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)), 
     edge.color = edge.cols.fall$col, add = T)

plot(fall.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
     vertex.shape = meta.fall.ab$shape_single, vertex.color = meta.fall.ab$shape_colour_single,
     edge.arrow.width = 0,edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)), add = T)

plot(fall.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
     vertex.shape = meta.fall.ab$shape_multiple, vertex.pie = meta.fall.ab$num.reg.ab.vector,
     vertex.pie.color = reg.ab.palette,edge.width = 0,edge.arrow.size = 0, edge.arrow.width = 0,
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, vertex.label.dist = 30, add = T)

# legend("bottomleft", title = as.expression(bquote(bold("Breeding range region"))), legend = c("Western", "Central", "Eastern"),
#        col = c( "#0072B2","#D55E00", "#009E73"), #reg.ab.palette[[1]],
#        pch = 15, cex = 1, box.lwd = 0)

################################################################################
# Population proportion by site during the spring season 
################################################################################

# Create a dataframe with the proportion of individuals from each section of the blackpoll warbler's range in each node
spring.stat.ab <- spring.stat %>% filter(!(is.na(next.cluster))) %>% merge(spring.breed.ab[,c("geo_id", "ab.unit")], by = "geo_id")

spring.ab.by.origin <- spring.stat.ab %>% group_by(cluster, Breeding_region_MC) %>%
  summarize(region.ab.units = sum(ab.unit), region.n = n_distinct(geo_id)) %>% ungroup() %>%
  complete(cluster, Breeding_region_MC, fill = list(region.ab.units = 0)) %>%
  complete(cluster, Breeding_region_MC, fill = list(region.n = 0))

# Convert data from wide to long
spring.ab.by.origin <- spring.ab.by.origin %>% pivot_wider(names_from = Breeding_region_MC, values_from = c(region.ab.units, region.n)) %>%
  rename(prop.ab.central = `region.ab.units_Central Region`,
         prop.ab.eastern = `region.ab.units_Eastern Region`,
         prop.ab.western = `region.ab.units_Western Region`,
         prop.ab.northwestern = `region.ab.units_Northwestern Region`,
         n.central = `region.n_Central Region`,
         n.eastern = `region.n_Eastern Region`,
         n.western = `region.n_Western Region`,
         n.northwestern = `region.n_Northwestern Region`)

# Merge abundance data with the node metadata
meta.spring.ab <- merge(meta.spring, spring.ab.by.origin, by.x = "vertex", by.y = "cluster") 

# Also add a column with the overall abundance at each site, and the number of individual at each site 
spring.stat.ab.per.site <- spring.stat.ab %>% group_by(cluster) %>% 
  summarize(r.abundance.at.cluster = sum(ab.unit))

meta.spring.ab <- merge(meta.spring.ab, spring.stat.ab.per.site, by.x = "vertex", by.y = "cluster") %>% 
  rowwise() %>%
  mutate(n.individuals.at.cluster = sum(n.central, n.eastern, n.western, n.northwestern))  

#create a column that can be converted to a numeric vector
meta.spring.ab <- transform(meta.spring.ab, num.reg.ab.vector = asplit(cbind(prop.ab.central, prop.ab.eastern, prop.ab.western, prop.ab.northwestern), 1))

# Plot of proportional node use during the spring migration 

# Create vector of vertex shapes 
meta.spring.ab <- meta.spring.ab %>% rowwise() %>%
  mutate(shape_single = length(which(c(prop.ab.central, 
                                       prop.ab.eastern, 
                                       prop.ab.western, 
                                       prop.ab.northwestern)==0))) %>%
  ungroup() %>%
  mutate(shape_single = ifelse(shape_single == 3 & node.type != "Breeding", "circle", "none")) %>%
  mutate(shape_single_breeding = ifelse(node.type == "Breeding", "square", "none"))%>%
  mutate(shape_multiple = ifelse(shape_single == "none" & shape_single_breeding == "none", "pie", "none")) %>%
  mutate(shape_colour_single = case_when(shape_single != "none" & prop.ab.central != 0 ~ "#0072B2",
                                         shape_single != "none" & prop.ab.eastern != 0 ~ "#009E73",
                                         shape_single != "none" & prop.ab.western != 0 ~ "#D55E00",
                                         shape_single != "none" & prop.ab.northwestern != 0 ~"#F0E442",
                                         .default = NA)) %>%
  mutate(shape_colour_single_breeding = case_when(shape_single_breeding != "none" & prop.ab.central != 0 ~ "#0072B2",
                                                  shape_single_breeding != "none" & prop.ab.eastern != 0 ~ "#009E73",
                                                  shape_single_breeding != "none" & prop.ab.western != 0 ~ "#D55E00",
                                                  shape_single_breeding != "none" & prop.ab.northwestern != 0 ~ "#F0E442",
                                                  .default = NA))

# Create a palette for site use by range region  
reg.ab.palette <- list(c("#0072B2", "#009E73", "#D55E00", "#F0E442"))

# Prepare edge colours for spring and spring edges (spring edges should not appear in this plot)
edge.cols.spring <- spring.con.ab %>% mutate(col = case_when(
  edge.type == "spring" ~ adjustcolor("darkgray", alpha.f = 0.9),
  edge.type == "spring" ~ NA))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

plot(spring.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
     vertex.shape = meta.spring.ab$shape_single_breeding, vertex.color = meta.spring.ab$shape_colour_single_breeding,
     edge.arrow.width = 0,edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)), 
     edge.color = edge.cols.spring$col, add = T)

plot(spring.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
     vertex.shape = meta.spring.ab$shape_single, vertex.color = meta.spring.ab$shape_colour_single,
     edge.arrow.width = 0,edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)), add = T)

plot(spring.graph.weighed.ab, vertex.size = 500, vertex.size2 = 200,
     vertex.shape = meta.spring.ab$shape_multiple, vertex.pie = meta.spring.ab$num.reg.ab.vector,
     vertex.pie.color = reg.ab.palette,edge.width = 0,edge.arrow.size = 0, edge.arrow.width = 0,
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, vertex.label.dist = 30, add = T)


################################################################################
# In the fall: calculate the proportion of individual using a node as a stopover
# Or as a a nonbreeding node
################################################################################

# calculate the ab.units of birds using nodes as breeding, nonbreeding and stopover areas
fall.use.per.node <- fall.stat.ab %>% group_by(cluster, site_type, geo_id) %>%
  summarize(ab.units = sum(ab.unit)) %>%
  group_by(cluster, site_type) %>%
  summarize(sum.use = sum(ab.units)) %>% ungroup() %>%
  complete(cluster, site_type, fill = list(sum.use = 0))

# Convert data from wide to long
fall.use.per.node <- fall.use.per.node  %>% pivot_wider(names_from = site_type, values_from = sum.use) %>%
  rename(use.breeding.ab = Breeding,
         use.stopover.ab = Stopover,
         use.nonbreeding.ab = Nonbreeding)

# Merge use data with the node metadata
meta.fall.ab <- merge(meta.fall.ab, fall.use.per.node , by.x = "vertex", by.y = "cluster")

#create a column that can be converted to a numeric vector
#meta.fall.ab <- transform(meta.fall.ab, use.ab.vector = asplit(cbind(use.breeding.ab,  use.stopover.ab, use.nonbreeding.ab), 1))

################################################################################
# In the spring: calculate the proportion of individual using a node as a stopover
# Or as a a nonbreeding node
################################################################################

# calculate the ab.units of birds using nodes as breeding, nonbreeding and stopover areas
spring.use.per.node <- spring.stat.ab %>% group_by(cluster, site_type, geo_id) %>%
  summarize(ab.units = sum(ab.unit)) %>%
  group_by(cluster, site_type) %>%
  summarize(sum.use = sum(ab.units)) %>% ungroup() %>%
  complete(cluster, site_type, fill = list(sum.use = 0))

# Convert data from wide to long
spring.use.per.node <- spring.use.per.node  %>% pivot_wider(names_from = site_type, values_from = sum.use) %>%
  rename(use.breeding.ab = Breeding,
         use.stopover.ab = Stopover,
         use.nonbreeding.ab = Nonbreeding)

# Merge use data with the node metadata
meta.spring.ab <- merge(meta.spring.ab, spring.use.per.node , by.x = "vertex", by.y = "cluster")

#create a column that can be converted to a numeric vector
#meta.spring.ab <- transform(meta.spring.ab, use.ab.vector = asplit(cbind(use.breeding.ab,  use.stopover.ab, use.nonbreeding.ab), 1))

################################################################################
### In the fall migration, stopover and nonbreeding specific node compositions 
################################################################################

# Composition of nonbreeding nodes at the end of the migration 
fall.nbr.ab <- fall.stat.ab %>% group_by(geo_id) %>% filter(sitenum == last(sitenum)) %>%
  group_by(cluster, Breeding_region_MC) %>% summarize(r.ab = sum(ab.unit)) %>%
  pivot_wider(names_from = Breeding_region_MC, values_from = c(r.ab), values_fill = 0) %>%
  mutate(tot.abundance = sum(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`),
         reg.no = ifelse(length(which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) == 0)) == 3, "single reg", "multi reg"),
         single_reg = ifelse(reg.no == "single reg", c("Eastern.Region", "Northwestern.Region", "Western.Region", "Central.Region")[which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) != 0)], NA)) %>%
  merge(meta.fall[,c("vertex", "Lon.50.", "Lat.50.")], by.x= "cluster", by.y = "vertex")

# Composition of nodes accoutnign only for use as stopover 
fall.stp.ab <- fall.stat.ab %>% group_by(geo_id) %>% filter(site_type != "Nonbreeding") %>%
  group_by(cluster, Breeding_region_MC) %>% summarize(r.ab = sum(ab.unit)) %>%
  pivot_wider(names_from = Breeding_region_MC, values_from = c(r.ab), values_fill = 0) %>%
  mutate(tot.abundance = sum(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`),
         reg.no = ifelse(length(which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) == 0)) == 3, "single reg", "multi reg"),
         single_reg = ifelse(reg.no == "single reg", c("Eastern.Region", "Northwestern.Region", "Western.Region", "Central.Region")[which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) != 0)], NA)) %>%
  merge(meta.fall[,c("vertex", "Lon.50.", "Lat.50.")], by.x= "cluster", by.y = "vertex")
  
#plot nonbreeding node use 
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-100, -40),ylim = c(-5, 15)) +
  geom_scatterpie(cols = c("Eastern Region","Northwestern Region", "Western Region", "Central Region"), colour = "black", data = fall.nbr.ab, mapping = aes(x = Lon.50., y = Lat.50., r = tot.abundance*10))

#plot stopover node use 
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-100, -40),ylim = c(-5, 15)) +
  geom_scatterpie(cols = c("Eastern Region","Northwestern Region", "Western Region", "Central Region"), colour = "black", data = fall.stp.ab, mapping = aes(x = Lon.50., y = Lat.50., r = tot.abundance*10))

################################################################################
### In the spring migration, stopover and nonbreeding specific node compositions 
################################################################################

# Composition of nonbreeding nodes at the start of the migration
spring.nbr.ab <- spring.stat.ab %>% group_by(geo_id) %>% filter(site_type == "Nonbreeding") %>%
  group_by(cluster, Breeding_region_MC) %>% summarize(r.ab = sum(ab.unit)) %>%
  pivot_wider(names_from = Breeding_region_MC, values_from = c(r.ab), values_fill = 0) %>%
  mutate(tot.abundance = sum(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`),
         reg.no = ifelse(length(which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) == 0)) == 3, "single reg", "multi reg"),
         single_reg = ifelse(reg.no == "single reg", c("Eastern.Region", "Northwestern.Region", "Western.Region", "Central.Region")[which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) != 0)], NA)) %>%
  merge(meta.spring[,c("vertex", "Lon.50.", "Lat.50.")], by.x= "cluster", by.y = "vertex") 

# Composition of nodes accounting only for use as stopover 
spring.stp.ab <- spring.stat.ab %>% group_by(geo_id) %>% filter(site_type != "Nonbreeding") %>%
  group_by(cluster, Breeding_region_MC) %>% summarize(r.ab = sum(ab.unit)) %>%
  pivot_wider(names_from = Breeding_region_MC, values_from = c(r.ab), values_fill = 0) %>%
  mutate(tot.abundance = sum(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`),
         reg.no = ifelse(length(which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) == 0)) == 3, "single reg", "multi reg"),
         single_reg = ifelse(reg.no == "single reg", c("Eastern.Region", "Northwestern.Region", "Western.Region", "Central.Region")[which(c(`Eastern Region`, `Northwestern Region`, `Western Region`, `Central Region`) != 0)], NA)) %>%
  merge(meta.spring[,c("vertex", "Lon.50.", "Lat.50.")], by.x= "cluster", by.y = "vertex")

#plot nonbreeding node use 
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-100, -40),ylim = c(-5, 15)) +
  geom_scatterpie(cols = c("Eastern Region","Northwestern Region", "Western Region", "Central Region"), data = spring.nbr.ab, mapping = aes(x = Lon.50., y = Lat.50., r = tot.abundance*10))

#plot stopover node use 
ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-100, -40),ylim = c(-5, 15)) +
  geom_scatterpie(cols = c("Eastern Region","Northwestern Region", "Western Region", "Central Region"), colour = "black", data = spring.stp.ab, mapping = aes(x = Lon.50., y = Lat.50., r = tot.abundance*10))

################################################################################
# Export data from the network construction
################################################################################

# Save elements necessary to build Fall network
write_csv(fall.stat, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.locations.csv")
write_csv(fall.move, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.all.locations.csv")
write_graph(fall.graph.weighed.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.graph.edge.list.txt",
            format = c("edgelist"))
write.csv(dplyr::select(meta.fall.ab, !num.reg.ab.vector), "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.csv")
write.csv(fall.con.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.edge.weights.csv")
write.csv(fall.stat, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stationary.data.csv")
write.csv(fall.edge.df.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.intra.cluster.movements.csv")
write.csv(fall.nbr.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.nbr.node.composition.csv")
write.csv(fall.stp.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.stp.node.composition.csv")
write.csv(fall.breed.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.abundance.per.bird.csv")

# Save elements necessary to build spring network
write_csv(spring.stat, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.stationary.locations.csv")
write_csv(spring.move, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.all.locations.csv")
write_graph(spring.graph.weighed.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.graph.edge.list.txt",
            format = c("edgelist"))
write.csv(dplyr::select(meta.spring.ab, !num.reg.ab.vector), "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.csv")
write.csv(spring.con.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.edge.weights.csv")
write.csv(spring.stat, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.stationary.data.csv")
write.csv(spring.edge.df.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.intra.cluster.movements.csv")
write.csv(spring.nbr.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.nbr.node.composition.csv")
write.csv(spring.stp.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.stp.node.composition.csv")
write.csv(spring.breed.ab, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.abundance.per.bird.csv")

##################################################################################
#  Nonbreeding range sub network for the  fall migration 
##################################################################################

# Add a column with inter-cluster connections to our location dataset
fall.stat.sub <- fall.stat %>% filter(Lat.50. < 13) %>% group_by(cluster) %>%
  filter(length(cluster) >1) %>% ungroup() %>%
  mutate(next.cluster = case_when(
  lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
  .default = NA)) 

# Create a dataframe with the edge info 
fall.edge.df.sub <- fall.stat.sub %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                            Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                            Lat.2.5., Lat.97.5., EndTime,
                                            sitenum, duration, period,study.site,
                                            Range_region, NB_count, period,
                                            site_type) %>% filter(!is.na(next.cluster))

# calculate node locations as the median locations of the points included within 
# Note, we use the 50% quantile of the posterior distribution 
fall.node.lon.sub <- fall.stat.sub %>% group_by(cluster) %>% 
  summarize(node.lon = median(Lon.50.))

fall.node.lat.sub <- fall.stat.sub%>% group_by(cluster) %>% 
  summarize(node.lat = median(Lat.50.))

# Create a layout with node locations 
meta.fall.sub <- data.frame("vertex" = sort(unique(fall.stat.sub$cluster)), 
                        #data.frame("vertex" = seq(min(fall.edge.df$cluster), max(fall.edge.df$cluster)), 
                        "Lon.50." = fall.node.lon.sub$node.lon,
                        "Lat.50." = fall.node.lat.sub$node.lat)

fall.location.sub <- as.matrix(meta.fall.sub[, c("Lon.50.", "Lat.50.")])

# Create the fall network
fall.graph.sub <- graph_from_data_frame(fall.edge.df.sub, directed = T, vertices = meta.fall.sub)

# Colour palette for site type
type.palette <- rainbow(3)  

# plot the fall network over North and South America
plot(wrld_simpl, xlim = c(-80, -54),ylim = c(-0, 13))
plot(fall.graph.sub, vertex.size = 200, vertex.size2 = 200,
     edge.width = 1, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = fall.location.sub, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), add = T)

# Add node weights 
fall.edge.df.ab.sub <- merge(fall.edge.df.sub, dplyr::select(fall.breed.ab, geo_id, ab.unit, br.region.prop.total.population, br.polygon))

# list of connections weighed by abundance unit
# We add an edge type if optional spring edges were added to distinguish fall and spring edges 
fall.con.ab.sub <- fall.edge.df.ab.sub %>% group_by(geo_id) %>%
  group_by(cluster, next.cluster) %>% 
  reframe(weight = sum(ab.unit), weight.n = n())

# Create a dataframe with the proportion of individuals from each section of the blackpoll warbler's range in each node
fall.stat.ab.sub <- fall.stat.sub %>% merge(fall.breed.ab[,c("geo_id", "ab.unit")], by = "geo_id") 

fall.ab.by.origin.sub <- fall.stat.ab.sub %>% group_by(cluster, Breeding_region_MC) %>%
  summarize(region.ab.units = sum(ab.unit), region.n = n_distinct(geo_id)) %>% ungroup() %>%
  complete(cluster, Breeding_region_MC, fill = list(region.ab.units = 0)) %>%
  complete(cluster, Breeding_region_MC, fill = list(region.n = 0)) %>%
  group_by(cluster) %>% mutate (ab.weight = sum(region.ab.units ), n.weight = sum(region.n))

# Convert data from wide to long
fall.ab.by.origin.sub <- fall.ab.by.origin.sub %>% pivot_wider(names_from = Breeding_region_MC, values_from = c(region.ab.units, region.n)) %>%
  rename(prop.ab.central = `region.ab.units_Central Region`,
         prop.ab.eastern = `region.ab.units_Eastern Region`,
         prop.ab.western = `region.ab.units_Western Region`,
         prop.ab.northwestern = `region.ab.units_Northwestern Region`,
         n.central = `region.n_Central Region`,
         n.eastern = `region.n_Eastern Region`,
         n.western = `region.n_Western Region`,
         n.northwestern = `region.n_Northwestern Region`)

# Merge abundance data with the node metadata
meta.fall.ab.sub <- merge(meta.fall.sub, fall.ab.by.origin.sub,  by.x = "vertex", by.y = "cluster") 

# Create weighed subgraph 
fall.graph.sub.weighed <- graph_from_data_frame(fall.con.ab.sub, directed = T, vertices = meta.fall.ab.sub)

save(fall.graph.sub.weighed, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.sub.graph.R")
write.csv(meta.fall.ab.sub, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Fall.node.metadata.sub.csv")

##################################################################################
#  Nonbreeding range sub network for the  spring migration 
##################################################################################

# Add a column with inter-cluster connections to our location dataset
spring.stat.sub <- spring.stat %>% filter(Lat.50. < 13) %>% group_by(cluster) %>%
  filter(length(cluster) >1) %>% ungroup() %>%
  mutate(next.cluster = case_when(
    lead(cluster) != cluster & lead(geo_id) == geo_id ~ lead(cluster),
    .default = NA)) 

# Create a dataframe with the edge info 
spring.edge.df.sub <- spring.stat.sub %>% dplyr::select(cluster, next.cluster, geo_id, StartTime,
                                                    Lon.50., Lon.2.5., Lon.97.5., Lat.50.,
                                                    Lat.2.5., Lat.97.5., EndTime,
                                                    sitenum, duration, period,study.site,
                                                    Range_region, NB_count, period,
                                                    site_type) %>% filter(!is.na(next.cluster))

# calculate node locations as the median locations of the points included within 
# Note, we use the 50% quantile of the posterior distribution 
spring.node.lon.sub <- spring.stat.sub %>% group_by(cluster) %>% 
  summarize(node.lon = median(Lon.50.))

spring.node.lat.sub <- spring.stat.sub%>% group_by(cluster) %>% 
  summarize(node.lat = median(Lat.50.))

# Create a layout with node locations 
meta.spring.sub <- data.frame("vertex" = sort(unique(spring.stat.sub$cluster)), 
                            #data.frame("vertex" = seq(min(spring.edge.df$cluster), max(spring.edge.df$cluster)), 
                            "Lon.50." = spring.node.lon.sub$node.lon,
                            "Lat.50." = spring.node.lat.sub$node.lat)

spring.location.sub <- as.matrix(meta.spring.sub[, c("Lon.50.", "Lat.50.")])

# Create the spring network
spring.graph.sub <- graph_from_data_frame(spring.edge.df.sub, directed = T, vertices = meta.spring.sub)

# Colour palette for site type
type.palette <- rainbow(3)  

# plot the spring network over North and South America
plot(wrld_simpl, xlim = c(-80, -54),ylim = c(-0, 13))
plot(spring.graph.sub, vertex.size = 200, vertex.size2 = 200,
     edge.width = 1, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = spring.location.sub, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), add = T)

# Add node weights 
spring.edge.df.ab.sub <- merge(spring.edge.df.sub, dplyr::select(spring.breed.ab, geo_id, ab.unit, br.region.prop.total.population, br.polygon))

# list of connections weighed by abundance unit
# We add an edge type if optional spring edges were added to distinguish spring and spring edges 
spring.con.ab.sub <- spring.edge.df.ab.sub %>% group_by(geo_id) %>%
  group_by(cluster, next.cluster) %>% 
  reframe(weight = sum(ab.unit), weight.n = n())

# Create a dataframe with the proportion of individuals from each section of the blackpoll warbler's range in each node
spring.stat.ab.sub <- spring.stat.sub %>% merge(spring.breed.ab[,c("geo_id", "ab.unit")], by = "geo_id") 

spring.ab.by.origin.sub <- spring.stat.ab.sub %>% group_by(cluster, Breeding_region_MC) %>%
  summarize(region.ab.units = sum(ab.unit), region.n = n_distinct(geo_id)) %>% ungroup() %>%
  complete(cluster, Breeding_region_MC, fill = list(region.ab.units = 0)) %>%
  complete(cluster, Breeding_region_MC, fill = list(region.n = 0)) %>%
  group_by(cluster) %>% mutate (ab.weight = sum(region.ab.units ), n.weight = sum(region.n))

# Convert data from wide to long
spring.ab.by.origin.sub <- spring.ab.by.origin.sub %>% pivot_wider(names_from = Breeding_region_MC, values_from = c(region.ab.units, region.n)) %>%
  rename(prop.ab.central = `region.ab.units_Central Region`,
         prop.ab.eastern = `region.ab.units_Eastern Region`,
         n.central = `region.n_Central Region`,
         n.eastern = `region.n_Eastern Region`)

# Merge abundance data with the node metadata
meta.spring.ab.sub <- merge(meta.spring.sub, spring.ab.by.origin.sub,  by.x = "vertex", by.y = "cluster") 

# Create weighed subgraph 
spring.graph.sub.weighed <- graph_from_data_frame(spring.con.ab.sub, directed = T, vertices = meta.spring.ab.sub)

save(spring.graph.sub.weighed, file = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.sub.graph.R")
write.csv(meta.spring.ab.sub, "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/Spring.node.metadata.sub.csv")
