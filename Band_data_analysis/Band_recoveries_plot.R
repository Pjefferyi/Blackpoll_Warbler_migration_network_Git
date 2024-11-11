
band.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Banding_data/NABBP_2022_grp_54.csv")

blpw.bands <- band.data %>% filter(SPECIES_ID == 6610)

#Plot blackpoll warbler band recoveries 
plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-10, 80))
points(blpw.bands$LON_DD, blpw.bands$LAT_DD)

blpw.bands.y <- blpw.bands %>% group_by(EVENT_YEAR) %>% summarize(encounters = n())

ggplot(data = blpw.bands.y, mapping = aes(x = EVENT_YEAR, y = encounters)) + 
  geom_point()
  