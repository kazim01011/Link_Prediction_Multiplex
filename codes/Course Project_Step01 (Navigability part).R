library(readxl)
library(muxViz)
library(igraph)
library(ggplot2)
library(openxlsx)
library(Matrix)
library(MASS)  
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(tictoc)
library(multinet)

# Define file paths
file_paths <- c(
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_1.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_2.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_3.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_4.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_5.xlsx"
)

# Read and store the data from each file
extended_edge_lists <- lapply(file_paths, function(fp) {
  network_data <- read_excel(fp)
  colnames(network_data) <- c("node.from", "layer.from", "node.to", "layer.to", "weight")
  network_data
})

# Combine the extended edge lists from all files
combined_extended_edge_list <- do.call(rbind, extended_edge_lists)

# Determine the number of layers and nodes
Layers <- length(file_paths)
Nodes  <- max(combined_extended_edge_list$node.from, combined_extended_edge_list$node.to)

edges <- transform(
  combined_extended_edge_list,
  from = node.from + Nodes * (layer.from - 1),
  to = node.to + Nodes * (layer.to - 1)
)

# Manually create the supra-adjacency matrix
SupraA <- sparseMatrix(
  i = edges$from,
  j = edges$to,
  x = edges$weight,
  dims = c(Nodes * Layers, Nodes * Layers)
)
SupraT <- BuildSupraTransitionMatrixFromSupraAdjacencyMatrix(SupraA, Layers, Nodes, Type="classical")

# Coverage calculation
TAUS <- 10^seq(-1,3,0.05)
cov.mul.appr <- GetCoverageEvolutionMultilayer(SupraT, Layers, Nodes, TAUS, Approximate = TRUE, Approximate.disconnected = igraph::clusters(igraph::graph.adjacency(SupraA))$no)
# Extract the real parts of the coverage values for the original network
cov_data <- transform(cov.mul.appr, Coverage = Re(rho), Network = "Original")

# Add new links and recalculate coverage
# Jaccard Links
jaccard_new_links <- list(
  c(176, 1, 174, 1, 108.595865),  # Node U, Layer U, Node V, Layer V, Weight
  c(2, 1, 36, 1, 1328.647309),
  c(2, 1, 11, 1, 514.877828),
  c(36, 1, 11, 1, 1602.568765),
  c(33, 5, 11, 5, 1054.989148),
  c(75, 5, 39, 5, 379.508026),
  c(133, 5, 141, 5, 30.926609)
)
#Adamic-Adar Links
adamic_adar_new_links <- list(
  c(176, 1, 174, 1, 108.595865),
  c(2, 1, 36, 1, 838.283119),
  c(2, 1, 11, 1, 324.851741),
  c(36, 1, 11, 1, 1011.108316),
  c(17, 5, 28, 5, 1108.454222),
  c(33, 5, 11, 5, 2109.978297),
  c(19, 5, 28, 5, 1512.438867),
  c(37, 5, 62, 5, 1933.825999),
  c(75, 5, 39, 5, 759.016053),
  c(133, 5, 141, 5, 30.926609),
  c(62, 5, 39, 5, 785.646357),
  c(62, 5, 20, 5, 986.263964),
  c(39, 5, 20, 5, 885.476380),
  c(11, 5, 28, 5, 1703.882906)
)

# Append new Jaccard links to Layer 1 and Layer 5
for (link in jaccard_new_links) {
  layer_index <- link[2]  # Layer index is the second element in each link vector
  extended_edge_lists[[layer_index]] <- rbind(extended_edge_lists[[layer_index]], data.frame(node.from = link[1], layer.from = link[2], node.to = link[3], layer.to = link[4], weight = link[5]))
}
for (link in adamic_adar_new_links) {
  layer_index <- link[2]  # Layer index is the second element in each link vector
  extended_edge_lists[[layer_index]] <- rbind(extended_edge_lists[[layer_index]], data.frame(node.from = link[1], layer.from = link[2], node.to = link[3], layer.to = link[4], weight = link[5]))
}

# Combine the updated extended edge lists
combined_extended_edge_list_updated <- do.call(rbind, extended_edge_lists)

# Rebuild the supra-adjacency matrix and recalculate coverage
edges_updated <- transform(
  combined_extended_edge_list_updated,
  from = node.from + Nodes * (layer.from - 1),
  to = node.to + Nodes * (layer.to - 1),
  weight = weight
)
SupraA2 <- sparseMatrix(
  i = edges_updated$from,
  j = edges_updated$to,
  x = edges_updated$weight,
  dims = c(Nodes * Layers, Nodes * Layers)
)
SupraT2 <- BuildSupraTransitionMatrixFromSupraAdjacencyMatrix(SupraA2, Layers, Nodes, Type="classical")
cov2.mul.appr <- GetCoverageEvolutionMultilayer(SupraT2, Layers, Nodes, TAUS, Approximate = TRUE, Approximate.disconnected = igraph::clusters(igraph::graph.adjacency(SupraA2))$no)

# Extract the real parts of the coverage values for the network with new links
cov_data2 <- transform(cov2.mul.appr, Coverage = Re(rho), Network = "With New Links")


# Ensure that both dataframes have the same columns and combine them
combined_cov_data <- rbind(cov_data, cov_data2)

# Plotting the comparison
ggplot(combined_cov_data, aes(x = tau, y = Coverage, color = Network)) +
  geom_line() +
  scale_x_log10() +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(
    title = "Coverage Evolution Comparison",
    x = "Diffusion Time (tau)",
    y = "Coverage",
    color = "Network"
  ) +
  theme(legend.position = "bottom")
