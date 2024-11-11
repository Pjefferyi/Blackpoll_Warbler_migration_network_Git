



# Time spent in th efall network 
fall.gplot.time.spent <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  #geom_text(data = fall.ggnet, mapping = aes(x = x, y = y, label = weight), nudge_x = 4)+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Time spent", 
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$time.spent, spring.ggnet$timespent)))+
  ggtitle("Fall network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(-0.18, 0.5), legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## Time spent spring network ----
spring.gplot.time.spent<- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, lwd = weight, colour = edge.type),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Time spent", 
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0),  max(fall.ggnet$time.spent, spring.ggnet$time.spent)))+
  ggtitle("Spring network") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## Second metric in fall network ----
fall.gplot.time.spent.ab <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_nodes(data = fall.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent.ab), shape=21, size  = 3)+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Relative use \nweighed by \ntime spent ", 
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0), max(fall.ggnet$time.spent.ab, spring.ggnet$time.spent.ab)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(-0.18, 0.5), legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))

## Second metric in spring network ----
spring.gplot.time.spent.ab  <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-170, -50),ylim = c(-3, 68)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  geom_nodes(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, fill = time.spent.ab), shape=21, size  = 3)+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  #geom_text(data = spring.ggnet, mapping = aes(x = x, y = y, cex = node.weight, label = participation.coef))+
  scale_fill_viridis_c(direction = -1, option = "magma", name = "Relative use \nweighed by \ntime spent", 
                       guide = guide_colorbar(frame.colour = "black"), limits = c(min(0), max(fall.gdata$time.spent.ab, fall.gdata$time.spent.ab)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "None", legend.key = element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        plot.title = element_text(size=10),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))


## create panel ----
metrics.time.fig <- (fall.gplot.time.spent | spring.gplot.time.spent)/ (fall.gplot.time.spent.ab| spring.gplot.time.spent.ab)

ggsave(plot = metrics.time.fig, filename = "nodes.metrics.time.png" ,  path = "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures", 
       units = "cm", width = 24*1.2, height = 10*1.2, dpi = "print", bg = "white")



# Figure 2: Node population composition ---- 

## fall population composition accounting for abundance during the fall migration  ----
fall.data <- igraph::as_data_frame(fall.graph, "vertices")
fall.gplot.comp.n <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = fall.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("black", alpha = 0.5), adjustcolor("blue", alpha = 0)), guide = "none")+
  geom_pie_glyph(slices= c("n.West", "n.Central", "n.East", "n.Northwest"), colour = "black", data = fall.data[fall.data$node.comp < 3,], mapping = aes(x = long, y = lat, radius = n.individuals)) +
  scale_radius(range = c(0.6, 3), unit = "mm", guide = "none")+
  scale_fill_manual(values = c("n.West" = "#D55E00", "n.Central" = "#009E73", "n.East" = "#0072B2", "n.Northwest" = "#F0E442"), name = "Breeding origin") +
  new_scale_fill()+
  geom_point(data = fall.data[fall.data$node.comp == 3,], mapping = aes(x = long, y = lat, fill = single.reg, size = n.individuals), shape= 21, colour = "black",  show.legend = F)+
  scale_size(range = c(1.2, 6), guide = "none")+
  scale_fill_manual(values = c("West" = "#D55E00", "Central" = "#009E73", "East" = "#0072B2", "Northwest" = "#F0E442"), name = "Breeding origin") +
  ggtitle("(a) Fall node use")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 10), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0,0,0,0), "pt"))+ 
  guides(fill = guide_legend(override.aes = list(size = 5)), )

## Spring population composition accounting for abundance during the spring migration ----
spring.data <- igraph::as_data_frame(spring.graph, "vertices")
spring.gplot.comp.n <- ggplot(st_as_sf(America))+
  geom_sf(colour = "black", fill = "#F7F7F7") +
  coord_sf(xlim = c(-165, -50),ylim = c(-5, 70)) +
  geom_edges(data = spring.ggnet, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = edge.type, lwd = weight),
             arrow = arrow(length = unit(9, "pt"), type = "closed", angle = 10))+
  scale_linewidth(range = c(0.1, 2), guide = "none")+
  scale_color_manual(values=c(adjustcolor("blue", alpha = 0), adjustcolor("black", alpha = 0.5)), guide = "none")+
  geom_pie_glyph(slices= c("n.West", "n.Central", "n.East", "n.Northwest"), colour = "black", data = spring.data, mapping = aes(x = long, y = lat, radius = n.individuals)) +
  scale_radius(range = c(0.6, 3), unit = "mm", guide = "none")+
  scale_fill_manual(values = c("n.West" = "#D55E00", "n.Central" = "#009E73", "n.East" = "#0072B2", "n.Northwest" = "#F0E442"), guide = "none") +
  new_scale_fill()+
  geom_point(data = spring.data[spring.data$node.comp == 3,], mapping = aes(x = long, y = lat, fill = single.reg, size = n.individuals), shape= 21, colour = "black")+
  geom_point(data = spring.data, mapping = aes(x = long, y = lat, fill = single.reg, size = node.weight), colour = NA, shape= 21)+
  scale_size(range = c(1.2, 7), breaks = c(0.1, 0.2, 0.3), name = "Node weight")+
  scale_fill_manual(values = c("West" = "#D55E00", "Central" = "#009E73", "East" = "#0072B2", "Northwest" = "#F0E442"), guide = "none") +
  ggtitle("(b) Spring node use")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.2, 0.4), text = element_text(size = 10), legend.key = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.line=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin= unit(c(0,0,0,0), "pt"))

node.comp.n.fig <- (fall.gplot.comp.n |spring.gplot.comp.n)


png("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Thesis_Documents/Figures/Node.comp.individuals.png", units = "cm", width = 25*1.2, height = 8*1.2, res = 500)
node.comp.n.fig
dev.off()