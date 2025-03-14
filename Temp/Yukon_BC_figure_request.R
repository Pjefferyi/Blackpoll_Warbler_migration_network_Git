# Plot the movements of blackpoll warblers from Yukon and British Columbia
library(sf)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(lubridate)
library(anytime)

# Load some polygons for countries, states, provinces, and the great lakes 
world_countries <- ne_countries(type = "countries", scale = "large")
world_states <- ne_states(country = c("United States of America", "Canada"))
America <- world_countries[((world_countries$continent %in%  c("North America", "South America") | world_countries$admin == "France") & world_countries$admin != "Greenland"),]
Lakes <- ne_download(scale = 50, type = "lakes", category = "physical")

# Load the dataset of bird movements
geo.all <- read.csv("C:/Users/Jelan/Github/Processed_blpw_movement_data/Temporary/BC.Yukon.locations.csv") %>% arrange(geo_id, StartTime) %>%
  group_by(geo_id) %>%
  mutate(period = ifelse(sitenum == 1 & geo_id != "WRMA04173", "Breeding", period),
         period = ifelse(sitenum == max(sitenum) & geo_id != "WRMA04173" & Recorded_North_South_mig == "South and partial North", "Failure", period),
         period = ifelse(sitenum == max(sitenum) & geo_id == "WRMA04173", "Breeding", period))

# Some birds have stopovers estimated close to their departure sites (within 250 km) - these are likely inaccurate, so we'll drop them
geo.all<- geo.all %>% group_by(geo_id) %>% filter(period != "Breeding" | sitenum == 1 | sitenum == max(sitenum))

# prepare the data to make it compatible with the  geom_segment() function 
# geom_segment allows us to highlight the portion of the geolocator location estimates for days within 2 weeks of the fall equinox 
geo.all.mod <- geo.all %>% mutate(next.Lat.50. = lead(Lat.50.), 
                                next.Lon.50. = lead(Lon.50.))%>%
  rowwise() %>%
  # calculate proximity to fall equinox
  mutate(StartTime = anytime(StartTime),
         EndTime = anytime(EndTime),
          equi_prox_start = difftime(StartTime, fall.equinox.date, units = "days"),
         equi_prox_end = difftime(EndTime, fall.equinox.date, units = "days"),
         two_weeks_equi_prox = ifelse(abs(equi_prox_start) <= 14 | abs(equi_prox_end) <= 14, "True", "False"),
  # calculate proximity to spring equinox (most birds migrate )
         equi_prox_spring_start = difftime(anytime(StartTime), spring.equinox.date, units = "days"),
         equi_prox_spring_end = difftime(EndTime, spring.equinox.date, units = "days"),
         two_weeks_equi_prox_spring = ifelse(abs(equi_prox_spring_start) <= 14 | abs(equi_prox_spring_end) <= 14, "True", "False")) 

# filter movement data for the fall migration of birds from Yukon and BC ----
ybc.fall.data <- geo.all.mod %>% filter(study.site %in% c("British Columbia, Highway 37", "Whitehorse, Yukon")) %>%
  group_by(geo_id) %>%
  filter(StartTime <= StartTime[which(NB_count == 1)])

# Here is the plot for the fall migration - you can uncomment the lines for the error bars
# This will show the 2.5 and 97.5% limits of the posterior distribution 
# but this will make the figure crowded
plot1 <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = world_states, aes(), colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-160, -40), ylim = c(-5, 65)) +
  # geom_errorbar(data = ybc.fall.data[ybc.fall.data$duration >=2,], aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  # geom_errorbar(data = ybc.fall.data[ybc.fall.data$duration >=2,], aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_segment(data = ybc.fall.data, mapping = aes(x = Lon.50., y = Lat.50., xend = next.Lon.50., yend = next.Lat.50., linetype = two_weeks_equi_prox, group = geo_id), col = "firebrick")+
  scale_linetype_manual(values = c("True" = 6, "False" = 1), labels = c("Migration route", "Low accuracy (equinox)"), name = NULL)+
  geom_point(data =  ybc.fall.data[ybc.fall.data$period %in% c("Breeding", "Non-breeding period"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period), cex = 3.5)+
  
  # Showing the location and duration spent at stopovers in this plot can be confusing due to the equinox
  # You can comment the two lines below to not do that
  geom_point(data =  ybc.fall.data[ybc.fall.data$period %in% c("Post-breeding migration") & ybc.fall.data$duration >=2,], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period, cex = duration))+
  scale_size("Stopover duration (days)")+
  
  scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4),
                     labels = c("Breeding site", "Nonbreeding site", "Post-breeding migration stopover", "Pre-breeding migration stopover"), name = NULL) +
  
  scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black"),
                      labels = c("Breeding site", "Nonbreeding site", "Post-breeding migration stopover", "Pre-breeding migration stopover"), name = NULL) +
  
  scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"),
                    labels = c("Breeding site", "Nonbreeding site", "Post-breeding migration stopover", "Pre-breeding migration stopover"), name = NULL)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        text = element_text(size = 12),
        legend.text=element_text(size=9),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        legend.position = c(0.18, 0.3), 
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.box.background= element_rect(fill = "white", colour = "black"))+
  guides(colour = guide_legend(order=1),
         fill = guide_legend(order=1),
         shape = guide_legend(order=1),
         linetype = guide_legend(order=2),
         size = guide_legend(order=3))

# Save the figure to a file path on your device
ggsave(plot = plot1, filename = "Western_Canada_fall_migration.jpg" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures",
       units = "cm", width = 15*1.2, height = 13*1.2, dpi = 500, bg = "white")

# filter movement data for the spring migration of birds from Yukon and BC ----
west.spring.data <- geo.all.mod %>% filter(study.site %in% c("British Columbia, Highway 37", "Whitehorse, Yukon")) %>%
  group_by(geo_id) %>%
  filter(StartTime >= StartTime[which(NB_count == max(NB_count, na.rm = T))] & Recorded_North_South_mig %in% c("Both"),)

# Here is the plot for the spring migration - you can again uncomment the lines for the error bars
plot2 <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  geom_sf(data = world_states, aes(), colour = "black", fill = "#F7F7F7") +
  geom_sf(data = Lakes, fill = "lightblue", lwd = 0.2, alpha = 1) +
  coord_sf(xlim = c(-160, -40), ylim = c(-5, 65)) +
  #geom_errorbar(data = west.spring.data, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  #geom_errorbar(data = west.spring.data, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.5, width = 1, alpha = 0.8, color = "black") +
  geom_segment(data = west.spring.data, mapping = aes(x = Lon.50., y = Lat.50., xend = next.Lon.50., yend = next.Lat.50., group = geo_id, linetype = two_weeks_equi_prox_spring,), col = "firebrick")+
  scale_linetype_manual(values = c("True" = 6, "False" = 1), breaks = c('False'), labels = c("Migration route"), name = NULL)+
  geom_point(data =  west.spring.data[west.spring.data$sitenum >= 2 & west.spring.data$period %in% c("Pre-breeding migration"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period,  cex = duration))+
  geom_point(data =  west.spring.data[west.spring.data$period %in% c("Breeding", "Non-breeding period"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = period, pch = period, col = period), cex = 3.5)+
  scale_shape_manual(values=c("Post-breeding migration" = 21 , "Non-breeding period"  = 22, "Pre-breeding migration" = 21, "Breeding"  = 24, "Failure" = 4),
                     labels = c("Breeding site", "Nonbreeding site", "Pre-breeding migration stopover", "post-breeding migration stopover"), name = NULL) +
  scale_colour_manual(values=c("Post-breeding migration" = "black" , "Non-breeding period"  = "white", "Pre-breeding migration" = "black", "Breeding"  = "white", "Failure" = "black"),
                      labels = c("Breeding site", "Nonbreeding site", "Pre-breeding migration stopover", "post-breeding migration stopover"), name = NULL) +
  scale_fill_manual(values=c("Post-breeding migration" = "#FDE725FF" , "Non-breeding period"  = "black", "Pre-breeding migration" = "#21908CFF", "Breeding"  = "black"),
                    labels = c("Breeding site", "Nonbreeding site", "Pre-breeding migration stopover", "post-breeding migration stopover"), name = NULL)+
  scale_size("Stopover duration (days)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        text = element_text(size = 12),
        legend.text=element_text(size=9),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.spacing = unit(-5, "pt"),
        legend.position = c(0.18, 0.28), 
        plot.margin = unit(c(0,0,0,0), "pt"),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.box.background= element_rect(fill = "white", colour = "black"))+
  guides(colour = guide_legend(order=1),
         fill = guide_legend(order=1),
         shape = guide_legend(order=1),
         linetype = guide_legend(order=2),
         size = guide_legend(order=3))

# Save the figure to a file path on your device
ggsave(plot = plot2, filename = "Western_Canada_spring_migration.jpg" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures",
       units = "cm", width = 15*1.2, height = 13*1.2, dpi = 500, bg = "white")
