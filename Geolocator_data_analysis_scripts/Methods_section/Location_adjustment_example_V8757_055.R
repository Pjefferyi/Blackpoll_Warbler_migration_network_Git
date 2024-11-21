# Plot the data from geolocator V8757_055 (Northwest territories) for an example of location adjustment for two stopovers at the same longitude but different latitudes 

dat <- findLocData(geo.ids = c("V8757_055"), check_col_length = F)

dat <- dat[(dat$sitenum > 0),]


range_crop <- st_crop(spData::world[(spData::world$continent %in% c("North America", "South America")),],
                      xmin = -140, xmax = -45, ymin = -20, ymax = 80)

no_adj_plot <- ggplot(range_crop) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = dat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") +
  geom_errorbar(data = dat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") +
  geom_point(data = dat, mapping = aes(x = Lon.50., y = Lat.50., size = duration), color = "blue", ) +
  geom_path(data = dat, mapping = aes(x = Lon.50., y = Lat.50.), color = "blue",  linewidth = 0.4) +
  geom_text(data = dat, mapping = aes(x = Lon.50., y = Lat.50., label = round(sitenum, digits = 0)),
            color = "black", nudge_x = 5, nudge_y = 3, check_overlap = T) +
  theme_bw() +
  labs(size = "Number of days site \nwas occupied") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dat <- dat[(dat$sitenum > 0 & dat$sitenum != 5),]
dat$duration[5] <- dat$duration[5]  + 6.5273148
dat$sitenum <- seq(1,11)


adj_plot <- ggplot(range_crop) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = dat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") +
  geom_errorbar(data = dat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") +
  geom_point(data = dat, mapping = aes(x = Lon.50., y = Lat.50., size = duration), color = "blue", ) +
  geom_path(data = dat, mapping = aes(x = Lon.50., y = Lat.50.), color = "blue",  linewidth = 0.4) +
  geom_text(data = dat, mapping = aes(x = Lon.50., y = Lat.50., label = round(sitenum, digits = 0)),
            color = "black", nudge_x = 5, nudge_y = 3, check_overlap = T) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  labs(size = "Number of days site \nwas occupied") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# Show the convergence of latitude after the equinox ###########################
  
#load the raw initial path x0_r
load(file = paste0(dir,"/", geo.id, "_initial_path_raw.csv"))

#plot of light transitions throughout the fall migration 
par(mfrow=c(2,1))
plot(twl$Twilight, x0_r[,1], type = "o", ylab = "Longitude", xlab = "Time", cex.lab=1.05, cex.axis=1.05)
rect(anytime("2012-10-03 10:23:17"), min(x0_r[,1])-2, anytime("2012-10-09 23:02:37"), max(x0_r[,1])+2, col = alpha("blue", 0.4), lty=0)
rect(anytime("2012-10-11 10:21:03"), min(x0_r[,1])-2, anytime("2012-11-04 10:47:57"), max(x0_r[,1])+2, col = alpha("red", 0.4), lty=0)
plot(twl$Twilight, x0_r[,2], type = "o", ylab = "Latitude", xlab = "Time", cex.lab=1.05, cex.axis=1.05)
rect(anytime("2012-10-03 10:23:17"), min(x0_r[,2])-2, anytime("2012-10-09 23:02:37"), max(x0_r[,2])+2, col = alpha("blue", 0.4), lty=0)
rect(anytime("2012-10-11 10:21:03"), min(x0_r[,2])-2, anytime("2012-11-04 10:47:57"), max(x0_r[,2])+2, col = alpha("red", 0.4), lty=0)


ggplot() + 
  geom_path(aes(x = twl$Twilight, y = x0_r[,1])) +
  geom_point(aes(x = twl$Twilight, y = x0_r[,1])) +
  theme_bw() +
  xlab("Time ") + 
  ylab("Latitude") +
  labs(size = "Number of days site \nwas occupied") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# Plot the data from geolocator V8757_096 (Newfoundland) for an example of location adjustment for a point at an unlikely latitude 
dat <- findLocData(geo.ids = c("V8757_134"), check_col_length = F)
dat <- dat[(dat$sitenum > 0),]

range_crop <- st_crop(spData::world[(spData::world$continent %in% c("North America", "South America")),],
                      xmin = -90, xmax = -45, ymin = -5, ymax = 60)

no_adj_plot <- ggplot(range_crop) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = dat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") +
  geom_errorbar(data = dat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") +
  geom_point(data = dat, mapping = aes(x = Lon.50., y = Lat.50., size = duration), color = "blue", ) +
  geom_path(data = dat, mapping = aes(x = Lon.50., y = Lat.50.), color = "blue",  linewidth = 0.4) +
  geom_text(data = dat, mapping = aes(x = Lon.50., y = Lat.50., label = round(sitenum, digits = 0)),
            color = "black", nudge_x = 5, nudge_y = 3, check_overlap = T) +
  theme_bw() +
  labs(size = "Number of days site \nwas occupied") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dat[(dat$sitenum == 3),]$Lat.50. <- 43.875599
dat[(dat$sitenum == 3),]$Lat.2.5. <- 43.875599
dat[(dat$sitenum == 3),]$Lat.97.5. <- 43.875599

adj_plot <- ggplot(range_crop) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = dat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") +
  geom_errorbar(data = dat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") +
  geom_point(data = dat, mapping = aes(x = Lon.50., y = Lat.50., size = duration), color = "blue", ) +
  geom_path(data = dat, mapping = aes(x = Lon.50., y = Lat.50.), color = "blue",  linewidth = 0.4) +
  geom_text(data = dat, mapping = aes(x = Lon.50., y = Lat.50., label = round(sitenum, digits = 0)),
            color = "black", nudge_x = 5, nudge_y = 3, check_overlap = T) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  labs(size = "Number of days site \nwas occupied") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Show the convergence of latitude after the equinox ###########################

#load the raw initial path x0_r
load(file = paste0(dir,"/", "V8757_134", "_initial_path_raw.csv"))

#plot of light transitions throughout the fall migration 
par(mfrow=c(2,1), mar = c(4,4,4,4))
plot(twl$Twilight, x0_r[,1], type = "o", ylab = "Longitude", xlab = "Time", cex.lab=1.05, cex.axis=1.05)
rect(anytime("2012-10-12 10:12:12"), min(x0_r[,1])-2, anytime("2012-10-16 21:54:36"), max(x0_r[,1])+2, col = alpha("orange", 0.4), lty=0)
plot(twl$Twilight, x0_r[,2], type = "o", ylab = "Latitude", xlab = "Time", cex.lab=1.05, cex.axis=1.05)
rect(anytime("2012-10-12 10:12:12"), min(x0_r[,2])-2, anytime("2012-10-16 21:54:36"), max(x0_r[,2])+2, col = alpha("orange", 0.4), lty=0)


# Plot the data from geolocator V8757_096 (Newfoundland) for an example of location adjustment for a point at an unlikely latitude 
dat <- findLocData(geo.ids = c("V8757_134"), check_col_length = F)
dat <- dat[(dat$sitenum > 0),]

range_crop <- st_crop(spData::world[(spData::world$continent %in% c("North America", "South America")),],
                      xmin = -90, xmax = -45, ymin = -5, ymax = 60)

no_adj_plot <- ggplot(range_crop) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = dat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") +
  geom_errorbar(data = dat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") +
  geom_point(data = dat, mapping = aes(x = Lon.50., y = Lat.50., size = duration), color = "blue", ) +
  geom_path(data = dat, mapping = aes(x = Lon.50., y = Lat.50.), color = "blue",  linewidth = 0.4) +
  geom_text(data = dat, mapping = aes(x = Lon.50., y = Lat.50., label = round(sitenum, digits = 0)),
            color = "black", nudge_x = 5, nudge_y = 3, check_overlap = T) +
  theme_bw() +
  labs(size = "Number of days site \was occupied") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dat[(dat$sitenum == 3),]$Lat.50. <- 43.875599
dat[(dat$sitenum == 3),]$Lat.2.5. <- 43.875599
dat[(dat$sitenum == 3),]$Lat.97.5. <- 43.875599

adj_plot <- ggplot(range_crop) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = dat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") +
  geom_errorbar(data = dat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") +
  geom_point(data = dat, mapping = aes(x = Lon.50., y = Lat.50., size = duration), color = "blue", ) +
  geom_path(data = dat, mapping = aes(x = Lon.50., y = Lat.50.), color = "blue",  linewidth = 0.4) +
  geom_text(data = dat, mapping = aes(x = Lon.50., y = Lat.50., label = round(sitenum, digits = 0)),
            color = "black", nudge_x = 5, nudge_y = 3, check_overlap = T) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  labs(size = "Number of days site \nwas occupied") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Show the convergence of latitude after the equinox ###########################

# load the raw initial path x0_r
load(file = paste0(dir,"/", "V8757_134", "_initial_path_raw.csv"))

# plot of light transitions throughout the fall migration 
par(mfrow=c(2,1), mar = c(4,4,4,4))
plot(twl$Twilight, x0_r[,1], type = "o", ylab = "Longitude", xlab = "Time", cex.lab=1.05, cex.axis=1.05)
rect(anytime("2012-10-12 10:12:12"), min(x0_r[,1])-2, anytime("2012-10-16 21:54:36"), max(x0_r[,1])+2, col = alpha("orange", 0.4), lty=0)
plot(twl$Twilight, x0_r[,2], type = "o", ylab = "Latitude", xlab = "Time", cex.lab=1.05, cex.axis=1.05)
rect(anytime("2012-10-12 10:12:12"), min(x0_r[,2])-2, anytime("2012-10-16 21:54:36"), max(x0_r[,2])+2, col = alpha("orange", 0.4), lty=0)

# Plot geolocators with uncertainty in location over the southern US ##########

geo.south <- findLocData(geo.ids = c("V8757_019", "V8757_029", "4068_014", "blpw14"), check_col_length = F)
plotLocVec(data = geo.south, stati_only = T, timing = c("Post-breeding migration", "Non-breeding period"))


