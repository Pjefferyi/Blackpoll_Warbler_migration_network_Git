dat <- findLocData(geo.ids = c("V8757_055"), check_col_length = F)

dat <- dat[(dat$sitenum > 0),]
# dat <- dat[(dat$sitenum > 0 & dat$sitenum != 5),]
# dat$duration[5] <- dat$duration[5]  + 6.5273148
# dat$sitenum <- seq(1,11)


range_crop <- st_crop(spData::world[(spData::world$continent %in% c("North America", "South America")),],
                      xmin = -140, xmax = -45, ymin = -20, ymax = 80)

ggplot(range_crop) +
  geom_sf() +
  coord_sf() +
  geom_errorbar(data = dat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), width=1, color = "firebrick") +
  geom_errorbar(data = dat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), width=1, color = "firebrick") +
  geom_point(data = dat, mapping = aes(x = Lon.50., y = Lat.50., size = duration), color = "blue", ) +
  geom_path(data = dat, mapping = aes(x = Lon.50., y = Lat.50.), color = "blue",  linewidth = 0.6) +
  geom_text(data = dat, mapping = aes(x = Lon.50., y = Lat.50., label = round(sitenum, digits = 0)),
            color = "black", nudge_x = 5, nudge_y = 3, check_overlap = T) +
  theme_bw() +
  #labs(size = "Number of days site \nwas occupied")
  theme(legend.position = "None",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
