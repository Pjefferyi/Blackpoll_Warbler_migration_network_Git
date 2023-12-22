# Plot of of nonbreeding movements 

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
library(fmdates)

#source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

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
                                   "D"), check_col_length = F)



# Define stationary periods as breeding, stopover, or non-breeding ##########################

geo.all <- geo.all %>% group_by(geo_id) %>% mutate(site_type = case_when(
  (sitenum == 1 | sitenum == max(sitenum)) & Recorded_North_South_mig == "Both" ~ "Breeding",
  sitenum == 1 & Recorded_North_South_mig %in% c("South and partial North", "South" ) ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "North" ~ "Breeding",
  (sitenum == max(sitenum)) & Recorded_North_South_mig == "South and partial North" & path.elongated == T ~ "Breeding",
  period == "Non-breeding period" & (duration >= 14 | sitenum == 1 | sitenum == max(sitenum)) ~ "Nonbreeding",
  .default = "Stopover"))

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

# data preparation ####################################################

# extract Nonbreeding stationary sites
NB.stat <- geo.all %>% dplyr::filter(sitenum > 0, site_type %in% c("Nonbreeding"),
                   period %in% c("Non-breeding period"),
                   Recorded_North_South_mig %in% c("Both", "South and partial North")) %>%
  ungroup()

# Add a new column specifying the number of nonbreeding movements performed by each individuals and whether or not each individual performed nonbreeding movements
NB.stat <- NB.stat %>% group_by(geo_id) %>% mutate(no.nbrs = n()-1,
                                                   nbr.mover = ifelse(no.nbrs  > 0, "mover", "nonmover")) %>% ungroup()

# Modify the data so we have two points for each movement between two sites (which means that we will have duplicate points)
for (i in seq(2, nrow(NB.stat)-1)){
  if (NB.stat[i-1,]$geo_id == NB.stat[i,]$geo_id & NB.stat[i+1,]$geo_id == NB.stat[i,]$geo_id){
    next.group <- NB.stat[i,]
    NB.stat <- add_row(NB.stat, next.group)
  }
}

# sort data by the geo.id 
NB.stat <- NB.stat %>% arrange(geo_id, StartTime)
#View(NB.stat[,c("geo_id", "nbr.move.num")])

# Add a new column to keep track of nonbreeding movements 
NB.stat$nbr.move.group <- NA
n <- 1

for (i in seq(1, nrow(NB.stat)-1)){
  if(NB.stat[i+1,]$geo_id == NB.stat[i,]$geo_id & NB.stat[i+1,]$StartTime != NB.stat[i,]$StartTime){
    NB.stat$nbr.move.group[i] <- n 
  }else{
    NB.stat$nbr.move.group[i] <- n 
    n <- n+1 
  }
}
NB.stat$nbr.move.group[i+1] <- n

# Check the departure time for each nonbreeding movement (EndTime of the first row in group) to determine whether the movement occured close to the equinox 
NB.stat <- NB.stat %>%
   mutate(spring.equinox.date = anytime(spring.equinox.date),
          equinox.nbr.move = ifelse(EndTime > anytime(spring.equinox.date) - days(14) & EndTime < anytime(spring.equinox.date) + days(14),
                                    "equinox affected", "equinox free"),
          timing.nbr.move = ifelse(month(anytime(EndTime)) > 10,
                                   "fall.nbr.movement", "winter.nbr.movements"))%>%
  mutate(equinox.nbr.move = ifelse(lag(nbr.move.group) == nbr.move.group, dplyr::lag(equinox.nbr.move, default = "a"), equinox.nbr.move),
         timing.nbr.move = ifelse(lag(nbr.move.group) == nbr.move.group, dplyr::lag(timing.nbr.move, default = "b"), timing.nbr.move))

NB.stat$equinox.nbr.move[1] <- "equinox free"
NB.stat$timing.nbr.move[1] <- "winter.nbr.movements"

#View(NB.stat[,c("nbr.move.group","equinox.nbr.move","timing.nbr.move", "StartTime")])

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "white", lwd = 0.2) +
  coord_sf(xlim = c(-90, -45),ylim = c(-10, 15)) +
  geom_point(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = as.factor(nbr.mover)), cex = 3, shape = 21, stroke = NA) +
  scale_fill_manual(values = c("nonmover" = "#009E73", "mover" = adjustcolor("gray", alpha = 0)), name = "Locations")+
  #geom_text(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., label = geo_id), cex = 2)+
  geom_path(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = as.factor(nbr.move.group), 
            linetype = equinox.nbr.move, col = timing.nbr.move),
            arrow = arrow(end = "last", type = "closed", length = unit(0.1, "inches")), lwd = 0.8) +
  scale_color_manual(values = c("#E66100", "#5D3A9B"), name = "Timing", label = c("October-December", "January-May"))+
  scale_linetype_manual(values = c("dashed", "solid"), name = "Spring equinox", label = c("within 14 days", "not within 14 days"))+
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.12, 0.4), legend.background = element_rect(fill = "transparent"))
  







