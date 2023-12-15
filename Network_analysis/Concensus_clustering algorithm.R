# Community clustering for network analysis 

# load libraries
library(igraph)
library(tidyr)

# Load helper functions 
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Load spring data for use as an example 
spring.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.graph.edge.list.txt", directed = TRUE)

# Load spring graph node metadata 
meta.spring.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.node.metadata.csv")

# Load spring graph edge weights
spring.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/spring.edge.weights.csv")

# Add weights to the spring graph, and convert the spring graph to an undirected graph 
E(spring.graph)$weight <- spring.con.ab$weight
undirected.spring.graph <- as.undirected(spring.graph, mode = "collapse",
                                       edge.attr.comb = "sum")

spring.communities.lab <- cluster_label_prop(undirected.spring.graph, weights = E(undirected.spring.graph)$weight)

############################### clustering parameters ##########################

algiter = 1000 #number of times the algorithm should be repeated 
thresh = 0.7 # threshold 

############################### Part 1 #########################################

# Function to run the algorithm on the network n times 
runAlg <- function(iterations, graph){
  
  comb.dt <- data.frame(comb = rep(NA, algiter ))
  
  for (i in seq(1:iterations)){
    
    comms <- cluster_label_prop(graph)
    mem <- comms$membership
    
    comb.dt$comb[i] <- list(mem)
  }
  
  # convert the output dataframe to a matrix
  comb <- matrix(unlist(comb.dt), ncol=length(comb.dt$comb[[1]]), byrow=TRUE)
  
  return(comb)
}

comb <- runAlg(iterations = algiter, graph = as.undirected(spring.graph))  

############################## Part 2 ##########################################
rows <- list()

# generate adjacency matrix by looping through node combinations and summing the number of times that they are in the same cluster 
for (i in 1:vcount(spring.graph)){
  
  eq <- 1:vcount(spring.graph)
  
  for (j in 1:vcount(spring.graph)){
  
    eq[j] <- sum(comb[1:algiter,i] == comb[1:algiter,j], na.rm=TRUE)
  }
  
  rows[[i]] <- eq/algiter
}
  
# adjacency matrix 
ad.mat<- matrix(unlist(rows), ncol=vcount(spring.graph), byrow=TRUE)

######################### part 3 ###############################################

#Set diag elements of matrix to 0 (to avoid loop edges)
diag(ad.mat) <- 0

# check for any disconnnected vertices (all weights < threshold)
# connect these vertices to their neighbour with the highest weight 
for (i in (seq(1,nrow(ad.mat)))){
  
  if (!(T %in% ad.mat[i,] >= thresh)){
    
    maxima <- which(ad.mat[i,] == max(ad.mat[,i]))
    
    ad.mat[i, maxima] <- thresh
    ad.mat[maxima, i] <- thresh
  }
}

# Set values below trheshold to 0 
ad.mat.ini <- ad.mat
ad.mat[ad.mat < thresh] <- 0

# transform to matrix into graph 
ad.graph <- graph_from_adjacency_matrix(ad.mat, weighted = T)
plot(ad.graph)

#cluster graph
comms <- cluster_label_prop(as.undirected(ad.graph))
plot(comms, ad.graph)

# plot communities on map
spring.comm.pal <- rainbow(length(seq(1, max(comms$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[comms$membership], add = T)

###################### part 4 ##################################################

# Until the number of partitions returned by the clustering algorithm is greater than 0, repeat step 2 and 3

iter.graph <- ad.graph
P <- 1 
iter <- 0 

while (P > 1){
  
  # run community detection algorithm on concensus network n times 
  comb <- runAlg(iterations = algiter, graph = as.undirected(iter.graph))  
  P <- nrow(unique(comb))
  
  rows <- list()
  # generate adjacency matrix by looping through node combinations and summing the number of times that they are in the same cluster 
  for (i in 1:vcount(iter.graph )){
    eq <- 1:vcount(iter.graph )
    for (j in 1:vcount(iter.graph )){
      eq[j] <- sum(comb[1:algiter,i] == comb[1:algiter,j], na.rm=TRUE)
    }
    rows[[i]] <- eq/algiter
  }
  # adjacency matrix 
  ad.mat <- matrix(unlist(rows), ncol=vcount(iter.graph ), byrow=TRUE)
  
  #Set diag elements of matrix to 0 (to avoid loop edges)
  diag(ad.mat) <- 0 
  
  # check for any disconnnected vertices (all weights < threshold)
  # connect these vertices to their neighbour with the highest weight 
  for (i in (seq(1,nrow(ad.mat)))){
    
    if (!(T %in% ad.mat[i,] >= thresh)){
      
      maxima <- which(ad.mat[i,] == max(ad.mat[,i]))
      
      ad.mat[i, maxima] <- thresh
      ad.mat[maxima, i] <- thresh
    }
  }
  
  # Set values below trheshold to 0 
  ad.mat[ad.mat < thresh] <- 0
  
  # Create a graph from the adjacency matrix
  iter.graph <- graph_from_adjacency_matrix(ad.mat, weighted = T)
  
  #keep track of iterations 
  iter <- iter + 1 
}

#plot final community structure
#cluster graph
comms.f <- cluster_label_prop(as.undirected(iter.graph))
plot(comms.f, iter.graph )

# plot communities on map
spring.comm.pal <- rainbow(length(seq(1, max(comms.f$membership))))

plot(wrld_simpl, xlim = c(-170, -30), ylim = c(-15, 70))
plot(spring.graph, vertex.label = NA, vertex.size = 200, vertex.size2 = 200,
     edge.width = spring.con.ab$weight*30, edge.arrow.size = 0, edge.arrow.width = 0,  
     layout = as.matrix(meta.spring.ab[, c("Lon.50.", "Lat.50.")]), rescale = F, asp = 0, xlim = c(-170, -30),
     ylim = c(-15, 70), edge.curved = rep(c(-0.05, 0.05), nrow(spring.con.ab)),
     vertex.color = spring.comm.pal[comms.f$membership], add = T)
  


