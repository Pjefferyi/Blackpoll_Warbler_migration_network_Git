#library(reticulate)
devtools::install_github("rstudio/reticulate", ref = "86ebb56")
library(igraph)

Sys.which("python")

#create a python environment with the required packages
#conda_create(envname = "blackpoll", packages = c("leidenalg", "igraph"))

#Choose enviroment 
use_condaenv("blackpoll")
py_config()

# Install required packages 
#conda_install(env = "blackpoll", packages = "leidenalg", "igraph", "numpy")
conda_install(env = "blackpoll", packages = "leidenalg", forge = T)

#View packages in environment and import modules
py_list_packages(envname = "blackpoll")
import("numpy")
import("igraph")
import("leidenalg")

# the leiden package 
#install.packages("leiden")
library(leiden)

# create sampel plot
adjacency_matrix <- rbind(cbind(matrix(round(rbinom(400, 1, 0.8)), 20, 20),
                                matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.1)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.8)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.2)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.1)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.9)), 20, 20)))


rownames(adjacency_matrix) <- 1:60
colnames(adjacency_matrix) <- 1:60
graph_object <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
graph_object

plot(graph_object, vertex.color = "grey75")

adjacency_matrix <- igraph::as_adjacency_matrix(graph_object)

partition <- leiden(adjacency_matrix)
