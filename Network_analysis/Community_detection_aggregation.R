
type.dt <- data.frame(type = rep(NA, 1000), score = rep(NA, 1000))
type.n <- 1

# dataframe with possible types 
for (i in seq(1:10000)){
  
   #comms <- cluster_label_prop(undirected.fall.graph, weights = E(undirected.fall.graph)$weight)
   comms <- cluster_label_prop(as.undirected(undirected.spring.graph))
   mem <- comms$membership
   
  if (list(mem) %in% type.dt$type){
    
    n <- match(list(mem), type.dt$type)
    type.dt$score[n] <- type.dt$score[n] + 1
    
    } else {
      
   type.dt$type[type.n] <- list(mem)
   type.dt$score[type.n] <- 1
   type.dt$n_clusters[type.n] <- length(unique(mem))
   type.dt$modularity[type.n] <- modularity(comms)
   type.n <- type.n +1
    }
}

type.dt <- type.dt %>% filter(!is.na(type)) %>% arrange(desc(score))

#find specific membership
membership <- unlist(type.dt$type[which(type.dt$score == max(type.dt$score))])
membership <- unlist(type.dt$type[3])

# # plot fall communities on map
# fall.comm.pal <- rainbow(length(seq(1, max(membership))))
# 
# plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
# plot(fall.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
#      edge.width = fall.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
#      layout = as.matrix(meta.fall.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
#      ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(fall.con.ab)),
#      vertex.color = fall.comm.pal[membership], add = T)

# plot spring communities on map
spring.comm.pal <- rainbow(length(seq(1, max(membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[membership], add = T)

################################################################################

#list of combinations
comb.dt <- data.frame(comb_id = rep(NA, 1000), comb = rep(NA, 1000))
n = 1

# dataframe with all combinations of memberships generated  
for (i in seq(1:1000)){
  
  comms <- cluster_label_prop(as.undirected(undirected.fall.graph))
  mem <- comms$membership
  
  comb.dt$comb[n] <- list(mem)
  comb.dt$comb_clusters[n] <- length(unique(mem))
  comb.dt$comb_modularity[n] <- modularity(comms)
  comb.dt$comb_id[n] <- n
  n = n+1 
}

comb.dt
hist(comb.dt$comb_clusters, breaks = 50)
mean(comb.dt$comb_modularity)

# We do the same with rewired graphs
comb.dt <- data.frame(comb_id = rep(NA, 1000), comb = rep(NA, 1000))
n = 1

for (i in seq(1:1000)){
  
  comms <- cluster_label_prop(as.undirected(rewire(undirected.spring.graph, with = keeping_degseq (niter = 100))))
  mem <- comms$membership
  
  comb.dt$comb[n] <- list(mem)
  comb.dt$comb_clusters[n] <- length(unique(mem))
  comb.dt$comb_modularity[n] <- modularity(comms)
  comb.dt$comb_id[n] <- n
  n = n+1 
}

comb.dt
hist(comb.dt$comb_clusters)
mean(comb.dt$comb_modularity)


