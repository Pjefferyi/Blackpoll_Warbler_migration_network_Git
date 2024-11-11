
#Load the helper functions (and associated packages)
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis/Geolocator_analysis_helper_functions.R")

#path to reference file
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"   

# Extract location data for spring migratory routes
geo.all <- findLocData(geo.ids = c(#"V8757_010",
                              #"V8296_004",
                              #"V8296_005",
                              #"V8296_006",
                              #"V8757_055",
                              #"V8757_018",
                              ##"V8757_021",
                              #"V8296_015",
                              #"V8296_017",
                              #"V8296_026",
                              ##"V8296_025",
                              #"V8296_007",
                              #"V8296_008",
                              #"V8757_019",
                              #"V8757_096",
                              #"V8757_134",
                              #"V8757_029",
                              #"V8757_078",
                              "blw09",
                              "blpw12",
                              "3254_001",
                              "4068_014",
                              "blpw14",
                              ##"3254_003",
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
                              "4210_004"#,
                              ##"4210_006",
                              ##"4210_010",
                              #"A",
                              #"B",
                              #"C",
                              #"D",
                              #"E",
                              #"WRMA04173"
                            ), check_col_length = F, ref_path = ref_path)
                            
                            #remove any NA values
geo.all <- geo.all[(!is.na(geo.all$Lon.50.) & !is.na(geo.all$"Lat.50.")),]

# Breeding period following the migration (only for complete tracks!)
geo.all <- geo.all %>% group_by(geo_id) %>%
  mutate(period = ifelse(sitenum == max(sitenum) & (Recorded_North_South_mig == "Both" | geo_id == "WRMA04173"), "Breeding period", period)) %>%
  ungroup()

# Prepare the data for the spring (only stationary locations) 
geo.spr.all <- geo.all[(geo.all$period %in% c("Pre-breeding migration", "Non-breeding period", "Breeding period")),]
geo.spr.stat <- geo.all[(geo.all$period %in% c("Pre-breeding migration", "Non-breeding period", "Breeding period") & geo.all$sitenum > 0),]

#Breeding and nonbreeding areas (not use for plotting) 
# Br.areas <- geo.all[(geo.all$sitenum == 1),]
# Br.areas[which(Br.areas$geo_id == "WRMA04173"), "Lon.50."] <- geo.all[(geo.all$geo_id == "WRMA04173"),]$Lon.50.[nrow(geo.all[(geo.all$geo_id == "WRMA04173"),])] 
# Br.areas[which(Br.areas$geo_id == "WRMA04173"), "Lat.50."] <- geo.all[(geo.all$geo_id == "WRMA04173"),]$Lat.50.[nrow(geo.all[(geo.all$geo_id == "WRMA04173"),])]
# Nbr.areas <- geo.all[(geo.all$period %in% c("Non-breeding period") & geo.all$duration > 50),]

#set the colours for each region
region.colours <- c("Churchill, Manitoba" = "pink", "Denali, Alaska" = "#D55E00", "Nome, Alaska" = "#009E73", "Whitehorse, Yukon" = "#0072B2")
symbols <- c(22, 23, 21)
symb.size <- c(4, 4, 4)

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-8, 70)) +
  #geom_errorbar(data = geo.spr.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5., color = Range_region), width=1, alpha = 0.2) + 
  #geom_errorbar(data = geo.spr.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5., color = Range_region), width=1, alpha = 0.6) + 
  geom_path(data = geo.spr.all, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, color = study.site), alpha = 0.5) +
  geom_point(data = geo.spr.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, fill = study.site, pch = period, size = period), colour = "black") +
  scale_fill_manual(values= region.colours, labels=c("Nome","Churchill", "Denali", "Whitehorse"),
                    name = "Breeding region")+
  scale_color_manual(values = region.colours, labels=c("Nome","Churchill", "Denali", "Whitehorse"), 
                     guide = guide_legend(override.aes = list(shape = 21, size = 2)),
                     name = "Breeding region")+
  scale_shape_manual(values = symbols, labels=c('Breeding site', 'Nonbreeding site', 'spring stopover'),
                     name = "Stationary areas", guide = guide_legend(override.aes = list(size = 4)))+
  scale_size_manual(values = symb.size, guide = "none", name = "Stationary areas")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.2, 0.4)) +
  xlab("Longitude") + 
  ylab("Latitude")+
  theme(text = element_text(size = 20))#+
#geom_text(data = geo.spr.stat, aes(x = Lon.50., y = Lat.50., group = geo_id),
#        label = geo.spr.stat$geo_id, size = 3)