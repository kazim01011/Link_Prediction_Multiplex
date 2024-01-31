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

# Define layer mapping
layer_mapping <- c('Layer_1' = 1, 'Layer_2' = 2, 'Layer_3' = 3, 'Layer_4' = 4, 'Layer_5' = 5)

#Read and combine new links
read_and_combine_new_links <- function(jaccard_file_path, adamic_adar_file_path, layer_mapping) {
  jaccard_links_df <- read.csv(jaccard_file_path, stringsAsFactors = FALSE)
  adamic_adar_links_df <- read.csv(adamic_adar_file_path, stringsAsFactors = FALSE)
  combined_df <- rbind(jaccard_links_df, adamic_adar_links_df)
  
  combined_links_list <- lapply(1:nrow(combined_df), function(idx) {
    list(
      Node_U = combined_df$Node_U[idx],
      Layer = layer_mapping[combined_df$Layer[idx]],
      Node_V = combined_df$Node_V[idx],
      weight = combined_df$weight[idx]
    )
  })
  return(combined_links_list)
}

#Add new links to network
add_links_to_network <- function(extended_edge_lists, new_links, layer_mapping) {
  for (link in new_links) {
    link_df <- data.frame(
      node.from = link$Node_U, 
      layer.from = layer_mapping[link$Layer], 
      node.to = link$Node_V, 
      layer.to = layer_mapping[link$Layer], 
      weight = link$weight
    )
    extended_edge_lists[[layer_mapping[link$Layer]]] <- rbind(extended_edge_lists[[layer_mapping[link$Layer]]], link_df)
  }
  return(extended_edge_lists)
}

#Analyse network
process_network <- function(extended_edge_lists, title_suffix, num_layers) {
    combined_extended_edge_list <- do.call(rbind, extended_edge_lists)
    Layers <- num_layers
    Nodes  <- max(max(combined_extended_edge_list$node.from), max(combined_extended_edge_list$node.to))
    
    # Preparing the edge list for BuildSupraAdjacencyMatrixFromExtendedEdgelist
    mEdges <- transform(
      combined_extended_edge_list,
      node.from = node.from,
      layer.from = layer.from,
      node.to = node.to,
      layer.to = layer.to,
      weight = weight
    )
    colnames(mEdges) <- c("node.from", "layer.from", "node.to", "layer.to", "weight")
    
    # Using BuildSupraAdjacencyMatrixFromExtendedEdgelist function
    SupraA <- BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges, Layers, Nodes, TRUE)
  
    SupraT <- BuildSupraTransitionMatrixFromSupraAdjacencyMatrix(SupraA, Layers, Nodes, Type="pagerank")
    TAUS <- 10^seq(-1,3,0.05)
    DisconnectedNodes <- igraph::clusters(igraph::graph.adjacency(SupraA))$no
    cov.mul.appr <- GetCoverageEvolutionMultilayer(SupraT, Layers, Nodes, TAUS, Approximate = TRUE, Approximate.disconnected = DisconnectedNodes)
    
    cov_data <- transform(cov.mul.appr, Coverage = Re(rho), Network = title_suffix)
    return(cov_data)
} 

# Function to generate a combined coverage plot of original and a specific stage
generate_combined_stage_plot <- function(original_cov_data, stage_cov_data, stage_title) {
  combined_data <- rbind(original_cov_data, transform(stage_cov_data, Network = stage_title))
  ggplot(combined_data, aes(x = tau, y = Coverage, color = Network)) +
    geom_line() +
    scale_x_log10() +
    theme_minimal() +
    labs(
      title = paste("Coverage Evolution - Original vs", stage_title),
      x = "Diffusion Time (tau)",
      y = "Coverage",
      color = "Stage"
    )
}

# Read and process the original data
original_file_paths <- c(
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_1.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_2.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_3.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_4.xlsx",
  "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Filtered_Case_5.xlsx"
)

extended_edge_lists_original <- list()
for (i in 1:length(original_file_paths)) {
  network_data <- read_excel(original_file_paths[i])
  colnames(network_data) <- c("node.from", "layer.from", "node.to", "layer.to", "weight")
  extended_edge_lists_original[[i]] <- network_data
}

# Analyze and visualize the original network
num_layers <- length(original_file_paths)
original_cov_data <- process_network(extended_edge_lists_original, "Original", num_layers)
#original_cov_plot <- generate_coverage_plot(original_cov_data, "Original Network Coverage")
#print(original_cov_plot)

############# Step 01 Addition of Indivisual Layer Predicted Links#############
# Define file paths for Jaccard and Adamic-Adar CSVs for individual layers
jaccard_individual_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_jaccard_individual.csv"
adamic_adar_individual_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_adamic_adar_individual.csv"

# Read and combine new links for individual layers
new_links_individual <- read_and_combine_new_links(jaccard_individual_path, adamic_adar_individual_path, layer_mapping)

# Add new links to the original network
extended_edge_lists_individual <- add_links_to_network(extended_edge_lists_original, new_links_individual, layer_mapping)

# Analyze and plot network after adding individual layer links
individual_cov_data <- process_network(extended_edge_lists_individual, "Individual Layer Links Added", num_layers)
individual_combined_plot <- generate_combined_stage_plot(original_cov_data, individual_cov_data, "Individual Layer Links Added")
print(individual_combined_plot)

############# Step 02 Addition of Two Layer Combination Predicted Links#############  

# Define file paths for two layer combination CSVs
jaccard_two_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_jaccard_two_layer.csv"
adamic_adar_two_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_adamic_adar_two_layer.csv"

# Read and combine new links for two layers
new_links_two_layers <- read_and_combine_new_links(jaccard_two_layers_path, adamic_adar_two_layers_path, layer_mapping)

# Add new links to the network with individual layer links
extended_edge_lists_two_layers <- add_links_to_network(extended_edge_lists_original, new_links_two_layers, layer_mapping)

# Analyze and plot network after adding two layer combination links
two_layer_cov_data <- process_network(extended_edge_lists_two_layers, "Two Layer Combination Links Added", num_layers)
two_layer_combined_plot <- generate_combined_stage_plot(original_cov_data, two_layer_cov_data, "Two Layer Combination Links Added")
print(two_layer_combined_plot)

############# Step 03 Addition of Three Layer Combination Predicted Links#############

# Define file paths for three layer combination CSVs
jaccard_three_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_jaccard_three_layer.csv"
adamic_adar_three_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_adamic_adar_three_layer.csv"

# Read and combine new links for three layers
new_links_three_layers <- read_and_combine_new_links(jaccard_three_layers_path, adamic_adar_three_layers_path, layer_mapping)

# Add new links to the network with two layer combination links
extended_edge_lists_three_layers <- add_links_to_network(extended_edge_lists_original, new_links_three_layers, layer_mapping)

# Analyze and plot network after adding three layer combination links
three_layer_cov_data <- process_network(extended_edge_lists_three_layers, "Three Layer Combination Links Added", num_layers)
three_layer_combined_plot <- generate_combined_stage_plot(original_cov_data, three_layer_cov_data, "Three Layer Combination Links Added")
print(three_layer_combined_plot)

############# Step 04 Addition of All Layer Combination Predicted Links#############

# Combine all new links
all_new_links <- c(new_links_individual, new_links_two_layers, new_links_three_layers)

# Add all new links to the original network
extended_edge_lists_all <- add_links_to_network(extended_edge_lists_original, all_new_links, layer_mapping)

# Analyze and plot network after adding all combination links
all_layer_cov_data <- process_network(extended_edge_lists_all, "All Layers Combination Links Added", num_layers)
all_layer_combined_plot <- generate_combined_stage_plot(original_cov_data, all_layer_cov_data, "All Layers Combination Links Added")
print(all_layer_combined_plot)

############# Step 05 Final Combined Plot#############

# Final combined plot of all stages
all_stages_combined_data <- rbind(
  transform(original_cov_data, Stage = "Original"),
  transform(individual_cov_data, Stage = "Individual"),
  transform(two_layer_cov_data, Stage = "Two Layer"),
  transform(three_layer_cov_data, Stage = "Three Layer"),
  transform(all_layer_cov_data, Stage = "All Layers")
)
final_combined_plot <- ggplot(all_stages_combined_data, aes(x = tau, y = Coverage, color = Stage)) +
  geom_line() +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "Combined Coverage Evolution Across All Stages",
    x = "Diffusion Time (tau)",
    y = "Coverage",
    color = "Stage"
  )
print(final_combined_plot)

############# Step 06 Final Combined Plot#############

generate_final_combined_plot <- function(original_cov_data, individual_cov_data, two_layer_cov_data, three_layer_cov_data) {
  all_stages_combined_data <- rbind(
    transform(original_cov_data, Stage = "Original"),
    transform(individual_cov_data, Stage = "Individual"),
    transform(two_layer_cov_data, Stage = "Two Layer"),
    transform(three_layer_cov_data, Stage = "Three Layer")
  )
  
  ggplot(all_stages_combined_data, aes(x = tau, y = Coverage, color = Stage, linetype = Stage)) +
    geom_line() +
    scale_x_log10() +
    scale_color_manual(values = c("black", "blue", "red", "green")) +
    scale_linetype_manual(values = c("Original" = "solid", "Individual" = "dashed", "Two Layer" = "dashed", "Three Layer" = "dotdash")) +  # Changed "Three Layer" to "dotdash"
    theme_minimal() +
    labs(
      title = "Combined Coverage Evolution Across All Stages",
      x = "Diffusion Time (tau)",
      y = "Coverage",
      color = "Stage",
      linetype = "Stage"
    )
}

# Generate and print the final combined plot
final_combined_plot <- generate_final_combined_plot(original_cov_data, individual_cov_data, two_layer_cov_data, three_layer_cov_data)
print(final_combined_plot)
