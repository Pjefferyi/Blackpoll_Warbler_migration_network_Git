
band.data <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Banding_data/NABBP_2022_grp_54.csv")

blpw.bands <- band.data %>% filter(SPECIES_ID == 6610, EVENT_TYPE == "E")

#Plot blackpoll warbler band recoveries 
plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-10, 80))
points(blpw.bands$LON_DD, blpw.bands$LAT_DD)
