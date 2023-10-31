# Modified consensus clustering function for use in robin's community analysis 

# consensusClusterMod #############################################################

# Input graph must be weighed and undirected 

concensusClusterMod <- function(graph, weights = NULL){
  
  ############################### Part 1 #########################################
  
  thresh <- 0.5 
  algiter <- 3000
  
  # Function to run the algorithm on the network n times 
  runAlg <- function(iterations, i.graph){
    
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
  
  comb <- runAlg(iterations = algiter, i.graph = as.undirected(graph))  
  
  ############################## Part 2 ##########################################
  rows <- list()
  
  # generate adjacency matrix by looping through node combinations and summing the number of times that they are in the same cluster 
  for (i in 1:vcount(graph)){
    
    eq <- 1:vcount(graph)
    
    for (j in 1:vcount(graph)){
      
      eq[j] <- sum(comb[1:algiter,i] == comb[1:algiter,j], na.rm=TRUE)
    }
    
    rows[[i]] <- eq/algiter
  }
  
  # adjacency matrix 
  ad.mat <- matrix(unlist(rows), ncol=vcount(graph), byrow=TRUE)
  
  ######################### part 3 ###############################################
  
  # Set values below trheshold to 0 
  ad.mat[ad.mat < thresh] <- 0
  
  #Set diag elements of matrix to 0 (to avoid loop edges)
  diag(ad.mat) <- 0
  
  # transform to matrix into graph 
  ad.graph <- graph_from_adjacency_matrix(ad.mat, weighted = T)
  
  ###################### part 4 ##################################################
  
  # Until the number of partitions returned by the clustering algorithm is greater than 0, repeat step 2 and 3
  
  iter.graph <- ad.graph
  P <- 1 
  iter <- 0 
  
  while (P > 1){
    
    # run community detection algorithm on consensus network n times 
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
    
    # Set values bewlo trheshold to 0 
    ad.mat[ad.mat < thresh] <- 0
    
    #Set diag elements of matrix to 0 (to avoid loop edges)
    diag(ad.mat) <- 0 
    
    iter.graph <- graph_from_adjacency_matrix(ad.mat, weighted = T)
    
    #keep track of iterations 
    iter <- iter + 1 
  }
  
  return(cluster_label_prop(iter.graph))
}

# Test calls for consensusClusterMod ##############################################

# Load fall data for use as an example
fall.graph <- read_graph("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/fall.graph.edge.list.txt", directed = TRUE)

# Load fall graph node metadata
meta.fall.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/fall.node.metadata.csv")

# Load fall graph edge weights
fall.con.ab <- read.csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Network_construction/fall.edge.weights.csv")

# Add weights to the fall graph, and convert the fall graph to an undirected graph
E(fall.graph)$weight <- fall.con.ab$weight
undirected.fall.graph <- as.undirected(fall.graph, mode = "collapse",
                                         edge.attr.comb = "sum")

# Run concensusCluster function
cluster_output <- concensusClusterMod(graph = undirected.fall.graph)
