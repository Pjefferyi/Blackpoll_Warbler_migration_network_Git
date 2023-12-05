# figures 

# Unweighed network sample --- 

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/unweighed_network_sample.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

# Plot the fall graph 
data("wrld_simpl")
type.palette <- viridis(3)  

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-80, -60), ylim = c(-4, 65), lwd = 0.5)

plot(fall.graph, vertex.size = 300, vertex.size2 = 200,
     edge.width = 0.5, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall.ab$node.type.num], vertex.label.dist = 30,
     add = T, vertex.label = NA, edge.color = "gray")
box(col = "black")

dev.off()

# Distribution in the Breeding grounds 

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/Breeding abundance.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

plot(bpw.fall.ab$breeding, xlim = c(-180, -50), ylim = c(35, 75), legend = F, axes = T, xlab = "Longitude", ylab = "Latitude")
plot(as_Spatial(abundance.regions), col = NA, border = "darkred", lwd = 2, add = T)
plot(wrld_simpl[wrld_simpl$NAME %in% c("United States", "Canada"),], add = T)
points(geo.breed$Lon.50., geo.breed$Lat.50., cex = 2, col = "blue", pch = 19)
#points(fall.stat$deploy.longitude, fall.stat$deploy.latitude, cex = 2, col = "blue", pch = 19)
#box(col = "black")

dev.off()

# Distribution in the Nonbreeding grounds

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/nonbreeding abundance.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

plot(bpw.fall.ab$nonbreeding, xlim = c(-100, -30), ylim = c(25, -15), legend = F, axes = T, xlab = "Longitude", ylab = "Latitude")
#plot(as_Spatial(abundance.regions), col = NA, border = "darkred", lwd = 2, add = T)
plot(wrld_simpl, add = T)
#points(geo.breed$Lon.50., geo.breed$Lat.50., cex = 2, col = "blue", pch = 19)
points( spring.stat[spring.stat$geo_id == "WRMA04173" & spring.stat$sitenum == 1,]$Lon.50., spring.stat[spring.stat$geo_id == "WRMA04173" & spring.stat$sitenum == 1,]$Lat.50., cex = 2, col = "blue", pch = 19)
#box(col = "black")

dev.off()


# Weighed network sample

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/weighed_network_sample.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

# Plot the fall graph 
data("wrld_simpl")
type.palette <- viridis(3)  

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-80, -60), ylim = c(-4, 65), lwd = 0.5)

plot(fall.graph, vertex.size = 300, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall.ab$node.type.num], vertex.label.dist = 30,
     add = T, vertex.label = NA, edge.color = adjustcolor("gray", alpha = 0.5))
box(col = "black")

dev.off()

# Stationary locations sample ----

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/clustering_sample.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = fall.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.2, alpha = 0.3, color = "black") +
  geom_errorbar(data = fall.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.2, alpha = 0.3, color = "black") +
  #geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5) +
  geom_point(data = fall.stat[(fall.stat$site_type!= "Breeding"),], mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 1) +
  labs(colour = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

dev.off()

# Sample of light data ---

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/light_data_sample.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

# We start with calibration based on the stationary periods before and after the migration
lightImage( tagdata = lig,
            offset = offset,     
            zlim = c(0, 20))

dev.off()

# Fall migratory network

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Travel/Colombia/Figures/Fall_migratory_network.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

# Plot the fall graph 
data("wrld_simpl")
type.palette <- type.palette <- c("#440154FF", "#FDE725FF", "#21908CFF")

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -20), ylim = c(-5, 60), col = "#F7F7F7", lwd = 0.5)

plot(as_Spatial(equinox_region), add = T, col = "#D9D5B2", lwd = 0.000000001)
#plot(as_Spatial(bpw_range), col = "#F4F0D3", lwd = 0.0001, add = T)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -20), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

plot(fall.graph, vertex.size = 400, vertex.size2 = 200,
     edge.width = fall.con.ab$weight*40, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = type.palette[meta.fall.ab$node.type.num], vertex.label.dist = 30,
     add = T, vertex.label = NA, edge.color = adjustcolor("darkgray", alpha.f = 0.6))

# legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
#        col = type.palette[unique(meta.fall.ab$node.type.num)],
#        pch = 19, title = "Node type",cex = 2)
# box(col = "black")

dev.off()

# Spring migratory network ---

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Travel/Colombia/Figures/spring_migratory_network.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

# Plot the spring graph 
data("wrld_simpl")
type.palette <- type.palette <- c("#440154FF", "#FDE725FF", "#21908CFF") 

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-165, -20), ylim = c(-5, 60), col = "#F7F7F7", lwd = 0.5)

plot(spring.graph, vertex.size = 400, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*40, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = type.palette[meta.spring.ab$node.type.num], vertex.label.dist = 30,
     add = T, vertex.label = NA, edge.color = adjustcolor("darkgray", alpha.f = 0.6))
# legend("bottomleft", legend = c("Stopover", "Nonbreeding", "Breeding"),
#        col = type.palette[unique(meta.spring.ab$node.type.num)],
#        pch = 19, title = "Node type",cex = 2)
# box(col = "black")

dev.off()

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/simple_network_legend.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

plot.new()
legend("bottomleft", legend = c("Breeding", "Stopover", "Nonbreeding"),
        col = c( "#440154FF", "#FDE725FF", "#21908CFF"),
       pch = 19, title = as.expression(bquote(bold("Node type"))),cex = 1, bty = "n")

dev.off()

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/node_population_composition_fall.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

# Create vector of vertex shapes 
meta.fall.ab <- meta.fall.ab %>% rowwise() %>%
  mutate(shape_single = length(which(c(prop.ab.central, 
                                       prop.ab.eastern, 
                                       prop.ab.western)==0))) %>%
  ungroup() %>%
  mutate(shape_single = ifelse(shape_single == 2 & node.type != "Breeding", "circle", "none")) %>%
  mutate(shape_single_breeding = ifelse(node.type == "Breeding", "square", "none"))%>%
  mutate(shape_multiple = ifelse(shape_single == "none" & shape_single_breeding == "none", "pie", "none")) %>%
  mutate(shape_colour_single = case_when(shape_single != "none" & prop.ab.central != 0 ~ "#009E73",
                                         shape_single != "none" & prop.ab.eastern != 0 ~ "#D55E00",
                                         shape_single != "none" & prop.ab.western != 0 ~ "#0072B2",
                                         .default = NA)) %>%
  mutate(shape_colour_single_breeding = case_when(shape_single_breeding != "none" & prop.ab.central != 0 ~ "#009E73",
                                                  shape_single_breeding != "none" & prop.ab.eastern != 0 ~ "#D55E00",
                                                  shape_single_breeding != "none" & prop.ab.western != 0 ~ "#0072B2",
                                                  .default = NA))

# Create a palette for site use by range region  
reg.ab.palette <- list(c("#009E73", "#D55E00", "#0072B2"))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

plot(fall.graph.weighed.ab, vertex.size = 300, vertex.size2 = 200,
     vertex.shape = meta.fall.ab$shape_single_breeding, vertex.color = meta.fall.ab$shape_colour_single_breeding,
     edge.arrow.width = 0,edge.width = fall.con.ab$weight *45, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

plot(fall.graph.weighed.ab, vertex.size = 450, vertex.size2 = 200,
     vertex.shape = meta.fall.ab$shape_single, vertex.color = meta.fall.ab$shape_colour_single,
     edge.arrow.width = 0,edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     edge.color = NA, add = T)

plot(fall.graph.weighed.ab, vertex.size = 450, vertex.size2 = 200,
     vertex.shape = meta.fall.ab$shape_multiple, vertex.pie = meta.fall.ab$num.reg.ab.vector,
     vertex.pie.color = reg.ab.palette,edge.width = 0,edge.arrow.size = 0, edge.arrow.width = 0,
     layout = fall.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label.dist = 30, vertex.label = NA,
     edge.color = NA, add = T)

dev.off()

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/node_population_composition_spring.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

# Create vector of vertex shapes 
meta.spring.ab <- meta.spring.ab %>% rowwise() %>%
  mutate(shape_single = length(which(c(prop.ab.central, 
                                       prop.ab.eastern, 
                                       prop.ab.west)==0))) %>%
  ungroup() %>%
  mutate(shape_single = ifelse(shape_single == 2 & node.type != "Breeding", "circle", "none")) %>%
  mutate(shape_single_breeding = ifelse(node.type == "Breeding", "square", "none"))%>%
  mutate(shape_multiple = ifelse(shape_single == "none" & shape_single_breeding == "none", "pie", "none")) %>%
  mutate(shape_colour_single = case_when(shape_single != "none" & prop.ab.central != 0 ~ "#009E73",
                                         shape_single != "none" & prop.ab.eastern != 0 ~ "#D55E00",
                                         shape_single != "none" & prop.ab.west != 0 ~ "#0072B2",
                                         .default = NA)) %>%
  mutate(shape_colour_single_breeding = case_when(shape_single_breeding != "none" & prop.ab.central != 0 ~ "#009E73",
                                                  shape_single_breeding != "none" & prop.ab.eastern != 0 ~ "#D55E00",
                                                  shape_single_breeding != "none" & prop.ab.west != 0 ~ "#0072B2",
                                                  .default = NA))

# Create a palette for site use by range region  
reg.ab.palette <- list(c("#009E73", "#D55E00", "#0072B2"))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = "#F7F7F7", lwd = 0.5)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -35), ylim = c(-10, 65), col = NA, lwd = 0.5, add = T)

plot(spring.graph.weighed.ab, vertex.size = 300, vertex.size2 = 200,
     vertex.shape = meta.spring.ab$shape_single_breeding, vertex.color = meta.spring.ab$shape_colour_single_breeding,
     edge.arrow.width = 0,edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

plot(spring.graph.weighed.ab, vertex.size = 450, vertex.size2 = 200,
     vertex.shape = meta.spring.ab$shape_single, vertex.color = meta.spring.ab$shape_colour_single,
     edge.arrow.width = 0,edge.width = 0, edge.arrow.size = 0, edge.arrow.width = 0,
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label = NA, edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     edge.color = NA, add = T)

plot(spring.graph.weighed.ab, vertex.size = 450, vertex.size2 = 200,
     vertex.shape = meta.spring.ab$shape_multiple, vertex.pie = meta.spring.ab$num.reg.ab.vector,
     vertex.pie.color = reg.ab.palette,edge.width = 0,edge.arrow.size = 0, edge.arrow.width = 0,
     layout = spring.location, rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), vertex.label.dist = 30, vertex.label = NA,
     edge.color = NA, add = T)

dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/node_composition_legend_color.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

plot.new()
legend("center", title = as.expression(bquote(bold("Breeding origin"))), legend = c("Western", "Central", "Eastern"),
       col = c("#0072B2", "#009E73", "#D55E00"),
       pch = 15, cex = 1, box.lwd = 0, bty = "n")

dev.off()

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/node_composition_legend_shape.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

plot.new()
legend("center", title = as.expression(bquote(bold("Node type"))), legend = c("Breeding", " Stopover or nonbreeding"),
       pch = c(0, 1), pt.cex= c(1.5, 2)*0.6, box.lwd = 0, cex = 0.5, bty = "n")

dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/spring_out_degree_centrality.jpg"),units = "cm", width = 12, height = 11,  quality = 100, res = 800)

## spring weighed degree centrality ----
spring.strength.c <- strength(spring.graph, mode = "all", weight = spring.con.ab$weight)

spring.con.ab.temp <- merge(spring.con.ab, meta.spring.ab[, c("vertex", "node.type")], by.x = "cluster", by.y = "vertex")

cust.color <- rep(adjustcolor("darkgray", alpha.f = 0.9), nrow(spring.con.ab.temp))
cust.color[spring.con.ab.temp$node.type != "Nonbreeding"] <- NA
cust.color[is.na(edge.cols.spring$col)]<- NA

cust.shape <- rep("circle", nrow(meta.spring.ab))
cust.shape[meta.spring.ab$node.type != "Nonbreeding"] <- "none"

# plot of the weighed centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(0, 1, 0.01)), palette = "Reds 3", rev = T) 
names(deg.c.palette) <- seq(0, 1, 0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-80, -38), ylim = c(-5, 15), col = "#F2F2F2", lwd = 0.5)

meta.spring.ab <- meta.spring.ab %>% mutate(breed.shape = ifelse(node.type == "Breeding", "square", "circle"))

plot(spring.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     vertex.shape = cust.shape,
     edge.width = spring.con.ab$weight*100, edge.arrow.size = 1, edge.arrow.width = 2,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = deg.c.palette[as.character(round(spring.strength.c, digits = 2))],
     edge.color = cust.color, add = T, vertex.label = NA)


dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/spring_in_degree_centrality.jpg"),units = "cm", width = 12, height = 11,  quality = 100, res = 800)

## spring weighed degree centrality ----
spring.strength.c <- strength(spring.graph, mode = "all", weight = spring.con.ab$weight)

spring.con.ab.temp <- merge(spring.con.ab, meta.spring.ab[, c("vertex", "node.type")], by.x = "next.cluster", by.y = "vertex") %>% 
  arrange(X)

cust.color <- rep(adjustcolor("darkgray", alpha.f = 0.9), nrow(spring.con.ab.temp))
cust.color[spring.con.ab.temp$node.type != "Nonbreeding"] <- NA

cust.shape <- rep("circle", nrow(meta.spring.ab))
cust.shape[meta.spring.ab$node.type != "Nonbreeding"] <- "none"

# plot of the weighed centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(0, 1, 0.01)), palette = "Reds 3", rev = T) 
names(deg.c.palette) <- seq(0, 1, 0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-80, -38), ylim = c(-5, 15), col = "#F2F2F2", lwd = 0.5)

meta.spring.ab <- meta.spring.ab %>% mutate(breed.shape = ifelse(node.type == "Breeding", "square", "circle"))

plot(spring.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     vertex.shape = cust.shape,
     edge.width = spring.con.ab$weight*50, edge.arrow.size = 0.3, edge.arrow.width = 3,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = deg.c.palette[as.character(round(spring.strength.c, digits = 2))],
     edge.color = cust.color, add = T, vertex.label = NA)


dev.off()



# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/fall_in_degree_centrality.jpg"),units = "cm", width = 12, height = 11,  quality = 100, res = 800)

## fall weighed degree centrality ----
fall.strength.c <- strength(fall.graph, mode = "all", weight = fall.con.ab$weight)

# plot of the weighed centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(0, 1, 0.01)), palette = "Reds 3", rev = T) 
names(deg.c.palette) <- seq(0, 1, 0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

meta.fall.ab <- meta.fall.ab %>% mutate(breed.shape = ifelse(node.type == "Breeding", "square", "circle"))

plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     vertex.shape = meta.fall.ab$breed.shape,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = deg.c.palette[as.character(round(fall.strength.c, digits = 2))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA)


dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/fall_out_degree_centrality.jpg"),units = "cm", width = 12, height = 11,  quality = 100, res = 800)

## fall weighed degree centrality ----
fall.strength.c <- strength(fall.graph, mode = "all", weight = fall.con.ab$weight)

fall.con.ab.temp <- merge(fall.con.ab, meta.fall.ab[, c("vertex", "node.type")], by.x = "cluster", by.y = "vertex")

cust.color <- rep(adjustcolor("darkgray", alpha.f = 0.9), nrow(fall.con.ab.temp))
cust.color[fall.con.ab.temp$node.type != "Nonbreeding"] <- NA
cust.color[is.na(edge.cols.fall$col)]<- NA

cust.shape <- rep("circle", nrow(meta.fall.ab))
cust.shape[meta.fall.ab$node.type != "Nonbreeding"] <- "none"

# plot of the weighed centrality of each node 
deg.c.palette <- hcl.colors(n = length(seq(0, 1, 0.01)), palette = "Reds 3", rev = T) 
names(deg.c.palette) <- seq(0, 1, 0.01)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-80, -38), ylim = c(-5, 15), col = "#F2F2F2", lwd = 0.5)

meta.fall.ab <- meta.fall.ab %>% mutate(breed.shape = ifelse(node.type == "Breeding", "square", "circle"))

plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     vertex.shape = cust.shape,
     edge.width = fall.con.ab$weight*100, edge.arrow.size = 1, edge.arrow.width = 2,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = deg.c.palette[as.character(round(fall.strength.c, digits = 2))],
     edge.color = cust.color, add = T, vertex.label = NA)

dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/degree centrality_legend.jpg"),units = "cm", width = 7.93, height = 8.21,  quality = 100, res = 800)

legend_image <- as.raster(matrix(hcl.colors(n = length(seq(0, 1, 0.01)), palette = "Reds 3", rev = F) , ncol=1))

plot(c(0,2),c(-0.01,1),type = 'n', axes = F,xlab = '', ylab = '')
title(main = 'Weighed degree centrality', cex.main = 0.9, line  = 0.8)
text(x= 1.25, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex = 1)
rasterImage(legend_image, 0.75, 0, 1,1)

dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/time_at_node_fall.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

time.palette <- hcl.colors(n = length(seq(0, max(meta.fall.ab$time.occupied +1))), palette = "Blues3", rev = T) 
names(time.palette) <- seq(0, max(meta.fall.ab$time.occupied +1))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = time.palette[as.character(round(meta.fall.ab$time.occupied))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA)

dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/time_at_node_spring.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

time.palette <- hcl.colors(n = length(seq(0, max(meta.spring.ab$time.occupied +1))), palette = "Blues3", rev = T) 
names(time.palette) <- seq(0, max(meta.spring.ab$time.occupied +1))

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),],
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)
plot(spring.graph, vertex.size = 300, vertex.size2 = 200, vertex.label.dist = 30,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = time.palette[as.character(round(meta.spring.ab$time.occupied))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T, vertex.label = NA)

dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/relative_abundance_node_spring.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

### estimated proportion of blackpoll warblers in each node 
r.ab.palette <- hcl.colors(n = length(seq(0, max(meta.spring.ab$r.abundance.at.cluster + 0.01), 0.001)), palette = "Reds 3", rev = T) 
names(r.ab.palette) <- seq(0, max(meta.spring.ab$r.abundance.at.cluster + 0.01), 0.001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(spring.graph, vertex.size = 300, vertex.size2 = 200, vertex.label= NA,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = r.ab.palette[as.character(round(meta.spring.ab$r.abundance.at.cluster, digits = 3))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

dev.off()

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/relative_abundance_node_fall.jpg"),units = "cm", width = 15.7, height = 20.1,  quality = 100, res = 800)

### estimated proportion of blackpoll warblers in each node
r.ab.palette <- hcl.colors(n = length(seq(0, max(meta.fall.ab$r.abundance.at.cluster + 0.01), 0.001)), palette = "Reds 3", rev = T) 
names(r.ab.palette) <- seq(0, max(meta.fall.ab$r.abundance.at.cluster + 0.01), 0.001)

plot(wrld_simpl[(wrld_simpl$REGION == 19 & wrld_simpl$NAME != "Greenland"),], 
     xlim = c(-170, -30), ylim = c(-15, 70), col = "#F2F2F2", lwd = 0.5)

plot(fall.graph, vertex.size = 300, vertex.size2 = 200, vertex.label= NA,
     edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
     vertex.color = r.ab.palette[as.character(round(meta.fall.ab$r.abundance.at.cluster, digits = 3))],
     edge.color = adjustcolor("darkgray", alpha.f = 0.6), add = T)

dev.off()


# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/Fall_clusters.jpg"),units = "cm", width = 15, height = 15,  quality = 100, res = 600)

fall.stat <- fall.stat %>% filter(site_type != "Breeding")

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = fall.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.3, alpha = 0.2, color = "black") +
  geom_errorbar(data = fall.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.3, alpha = 0.2, color = "black") +
  #geom_path(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5, linewidth = 0.2, linetype="dashed" ) +
  #geom_point(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 1.5, alpha = 0.8) +
  geom_point(data = fall.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), cex = 1.5, alpha = 0.2) +
  #geom_text(data = meta.fall[(meta.fall$node.type != "Breeding"), ], aes(x = Lon.50., y = Lat.50., label = vertex))+
  labs(colour = "Cluster", x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None")

dev.off()

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/spring_clusters.jpg"),units = "cm", width = 15, height = 15,  quality = 100, res = 600)

spring.stat <- spring.stat %>% filter(site_type != "Breeding")
#spring.stat.large <- spring.stat %>% group_by(cluster) %>% filter(n() > 3)

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = NA, fill = "lightgray") +
  coord_sf(xlim = c(-170, -30),ylim = c(-15, 70)) +
  geom_errorbar(data = spring.stat, aes(x = Lon.50., ymin= Lat.2.5., ymax= Lat.97.5.), linewidth = 0.3, alpha = 0.2, color = "black") +
  geom_errorbar(data = spring.stat, aes(y = Lat.50., xmin= Lon.2.5., xmax= Lon.97.5.), linewidth = 0.3, alpha = 0.2, color = "black") +
  #geom_path(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.5, linewidth = 0.2, linetype="dashed" ) +
  #geom_point(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id, colour = as.factor(cluster)), cex = 2, alpha = 0.8) +
  geom_point(data = spring.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), cex = 2, alpha = 0.2) +
  #geom_text(data = meta.spring[(meta.spring$node.type != "Breeding"), ], aes(x = Lon.50., y = Lat.50., label = vertex))+
  labs(colour = "Cluster", x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "None")
  #geom_density2d(data = spring.stat.large, aes(x = Lon.50., y = Lat.50., group = cluster, colour = as.factor(cluster)))

dev.off()

# open jpeg
jpeg(paste0("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Conferences/AOS_2023/AOS_poster/AOS_poster_figures/nonbreeding_movements.jpg"),units = "cm", width = 15, height = 15,  quality = 100, res = 600)

ggplot(st_as_sf(wrld_simpl))+
  geom_sf(colour = "black", fill = "#F7F7F7", lwd = 0.2) +
  coord_sf(xlim = c(-90, -45),ylim = c(-10, 15)) +
  geom_point(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), cex = 1.5, col = "firebrick") +
  geom_path(data = NB.stat, mapping = aes(x = Lon.50., y = Lat.50., group = geo_id), alpha = 0.7,
            arrow = arrow(end = "last", type = "open", length = unit(0.12, "inches")), lwd = 0.5, col = "blue") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dev.off()