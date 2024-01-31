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
process_network <- function(extended_edge_lists, title_suffix, num_layers,rw_type, directed) {
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
  SupraA <- BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges, Layers, Nodes, directed)
  
  SupraT <- BuildSupraTransitionMatrixFromSupraAdjacencyMatrix(SupraA, Layers, Nodes, Type = rw_type)
  TAUS <- 10^seq(-1,3,0.05)
  DisconnectedNodes <- igraph::clusters(igraph::graph.adjacency(SupraA))$no
  cov.mul.appr <- GetCoverageEvolutionMultilayer(SupraT, Layers, Nodes, TAUS, Approximate = TRUE, Approximate.disconnected = DisconnectedNodes)
  
  cov_data <- transform(cov.mul.appr, Coverage = Re(rho), Network = title_suffix)
  return(cov_data)
} 

# RW strategies and network types
rw_types <- c("pagerank", "diffusive", "maxent", "physical", "relaxed-physical")
directed_options <- c(TRUE, FALSE)

# Function to generate a combined coverage plot of original and a specific stage
generate_combined_stage_plot <- function(original_cov_data, stage_cov_data, stage_title) {
  combined_data <- rbind(original_cov_data, transform(stage_cov_data, Network = stage_title))
  ggplot(combined_data, aes(x = tau, y = Coverage, color = Network, linetype = Network)) +
    geom_line() +
    geom_point(aes(shape = Network), size = 1.5) +  # Adds markers with different shapes
    scale_x_log10() +
    scale_shape_manual(values = c(17, 19)) +  # 17: Triangle, 19: Circle
    scale_linetype_manual(values = c("solid", "dashed")) +
    theme_minimal() +
    labs(
      title = paste("Coverage Evolution - Original vs", stage_title),
      x = "Diffusion Time (tau)",
      y = "Coverage",
      color = "Network",
      linetype = "Network",
      shape = "Network"
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

num_layers <- length(original_file_paths)

# ---------------------------------------------------------------------------------------------------------------
                                        # Directed Network Analysis
# ---------------------------------------------------------------------------------------------------------------

# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, << --- PageRank --- >> ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Analyze and visualize the original network (Directed, PageRank)
original_cov_data_pagerank <- process_network(extended_edge_lists_original, "Original - Directed - PageRank", num_layers, "pagerank", TRUE)

# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Individual Layer Analysis ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

jaccard_individual_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_jaccard_individual.csv"
adamic_adar_individual_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_adamic_adar_individual.csv"

new_links_individual <- read_and_combine_new_links(jaccard_individual_path, adamic_adar_individual_path, layer_mapping)
extended_edge_lists_individual <- add_links_to_network(extended_edge_lists_original, new_links_individual, layer_mapping)

# Process the network with individual layer edges added (Directed, PageRank)
individual_cov_data <- process_network(extended_edge_lists_individual, "Individual Layer - Directed - PageRank", num_layers, "pagerank", TRUE)

# Generate and print the plot comparing the original network with the individual layer added network
individual_layer_plot <- generate_combined_stage_plot(original_cov_data_pagerank, individual_cov_data, "Individual Layer Links Added - PageRank")
print(individual_layer_plot)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Two Layer Addition,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Read and combine new links for two layer combinations
jaccard_two_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_jaccard_two_layer.csv"
adamic_adar_two_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_adamic_adar_two_layer.csv"

# Read and combine new links for two layers
new_links_two_layers <- read_and_combine_new_links(jaccard_two_layers_path, adamic_adar_two_layers_path, layer_mapping)

# Add new links to the network with individual layer links
extended_edge_lists_two_layers <- add_links_to_network(extended_edge_lists_original, new_links_two_layers, layer_mapping)

# Analyze and plot network after adding two layer combination links
two_layer_cov_data <- process_network(extended_edge_lists_two_layers, "Two Layer - Directed - PageRankd", num_layers, "pagerank", TRUE)
two_layer_combined_plot <- generate_combined_stage_plot(original_cov_data_pagerank, two_layer_cov_data, "Two Layer Combination Links Added")
print(two_layer_combined_plot)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Three Layer Addition,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Read and combine new links for three layer combinations

jaccard_three_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_jaccard_three_layer.csv"
adamic_adar_three_layers_path <- "/Volumes/Data/NDSU/PhD Work/Research/IME Research/AI-Energy/Data/Updated Data/Correct Predicted Data/cor_adamic_adar_three_layer.csv"

# Read and combine new links for three layers
new_links_three_layers <- read_and_combine_new_links(jaccard_three_layers_path, adamic_adar_three_layers_path, layer_mapping)

# Add new links to the network with two layer combination links
extended_edge_lists_three_layers <- add_links_to_network(extended_edge_lists_original, new_links_three_layers, layer_mapping)

# Analyze and plot network after adding three layer combination links
three_layer_cov_data <- process_network(extended_edge_lists_three_layers, "Three Layer Combination Links Added", num_layers, "pagerank", TRUE)
three_layer_combined_plot <- generate_combined_stage_plot(original_cov_data_pagerank, three_layer_cov_data, "Three Layer Combination Links Added")
print(three_layer_combined_plot)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, All Layer Addition,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Combine all new links
all_new_links <- c(new_links_individual, new_links_two_layers, new_links_three_layers)

# Add all new links to the original network
extended_edge_lists_all <- add_links_to_network(extended_edge_lists_original, all_new_links, layer_mapping)

# Analyze and plot network after adding all combination links
all_layer_cov_data <- process_network(extended_edge_lists_all, "All Layers Combination Links Added", num_layers, "pagerank", TRUE)
all_layer_combined_plot <- generate_combined_stage_plot(original_cov_data_pagerank, all_layer_cov_data, "All Layers Combination Links Added")
print(all_layer_combined_plot)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  Classical ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Analyze and visualize the original network (Directed, Classical)
original_cov_data_cl <- process_network(extended_edge_lists_original, "Original - Directed - Classical", num_layers, "classical", TRUE)

# Process the network with individual layer edges added (Directed, Classical)
individual_cov_data_cl <- process_network(extended_edge_lists_individual, "Individual Layer - Directed - Classical", num_layers, "classical", TRUE)

# Generate and print the plot comparing the original network with the individual layer added network
individual_layer_plot_cl <- generate_combined_stage_plot(original_cov_data_cl, individual_cov_data_cl, "Individual Layer Links Added - Classical")
print(individual_layer_plot_cl)

# Analyze and plot network after adding two layer combination links
two_layer_cov_data_cl <- process_network(extended_edge_lists_two_layers, "Two Layer - Directed - Classical", num_layers, "classical", TRUE)
two_layer_combined_plot_cl <- generate_combined_stage_plot(original_cov_data_cl, two_layer_cov_data_cl, "Two Layer Combination Links Added - Classical")
print(two_layer_combined_plot_cl)

# Analyze and plot network after adding three layer combination links
three_layer_cov_data_cl <- process_network(extended_edge_lists_three_layers, "Three Layer Combination Links Added - Classical", num_layers, "classical", TRUE)
three_layer_combined_plot_cl <- generate_combined_stage_plot(original_cov_data_cl, three_layer_cov_data_cl, "Three Layer Combination Links Added - Classical")
print(three_layer_combined_plot_cl)

# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Diffusive ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Analyze and visualize the original network (Directed, Diffusive)
original_cov_data_df <- process_network(extended_edge_lists_original, "Original - Directed - Diffusive", num_layers, "diffusive", TRUE)

# Process the network with individual layer edges added (Directed, Diffusive)
individual_cov_data_df <- process_network(extended_edge_lists_individual, "Individual Layer - Directed - Diffusive", num_layers, "diffusive", TRUE)

# Generate and print the plot comparing the original network with the individual layer added network
individual_layer_plot_df <- generate_combined_stage_plot(original_cov_data_df, individual_cov_data_df, "Individual Layer Links Added - Diffusive")
print(individual_layer_plot_df)

# Analyze and plot network after adding two layer combination links
two_layer_cov_data_df <- process_network(extended_edge_lists_two_layers, "Two Layer - Directed - Diffusive", num_layers, "diffusive", TRUE)
two_layer_combined_plot_df <- generate_combined_stage_plot(original_cov_data_df, two_layer_cov_data_df, "Two Layer Combination Links Added")
print(two_layer_combined_plot_df)

#Analyzing and plotting network after adding three layer combination links
three_layer_cov_data_df <- process_network(extended_edge_lists_three_layers, "Three Layer Combination Links Added", num_layers, "diffusive", TRUE)
three_layer_combined_plot_df <- generate_combined_stage_plot(original_cov_data_df, three_layer_cov_data_df, "Three Layer Combination Links Added")
print(three_layer_combined_plot_df)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Maxent ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Analyze and visualize the original network (Directed, Maxent)
original_cov_data_maxent <- process_network(extended_edge_lists_original, "Original - Directed - Maxent", num_layers, "maxent", TRUE)

# Process the network with individual layer edges added (Directed, Maxent)
individual_cov_data_maxent <- process_network(extended_edge_lists_individual, "Individual Layer - Directed - Maxent", num_layers, "maxent", TRUE)

# Generate and print the plot comparing the original network with the individual layer added network
individual_layer_plot_maxent <- generate_combined_stage_plot(original_cov_data_maxent, individual_cov_data_maxent, "Individual Layer Links Added - Maxent")
print(individual_layer_plot_maxent)

# Analyze and plot network after adding two layer combination links
two_layer_cov_data_maxent <- process_network(extended_edge_lists_two_layers, "Two Layer - Directed - Maxent", num_layers, "maxent", TRUE)
two_layer_combined_plot_maxent <- generate_combined_stage_plot(original_cov_data_maxent, two_layer_cov_data_maxent, "Two Layer Combination Links Added")
print(two_layer_combined_plot_maxent)

#Analyzing and plotting network after adding three layer combination links
three_layer_cov_data_maxent <- process_network(extended_edge_lists_three_layers, "Three Layer Combination Links Added", num_layers, "maxent", TRUE)
three_layer_combined_plot_maxent <- generate_combined_stage_plot(original_cov_data_maxent, three_layer_cov_data_maxent, "Three Layer Combination Links Added")
print(three_layer_combined_plot_maxent)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Physical ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Analyze and visualize the original network (Directed, Physical)
original_cov_data_physical <- process_network(extended_edge_lists_original, "Original - Directed - Physical", num_layers, "physical", FALSE)

# Process the network with individual layer edges added (Directed, Physical)
individual_cov_data_physical <- process_network(extended_edge_lists_individual, "Individual Layer - Directed - Physical", num_layers, "physical", TRUE)

# Generate and print the plot comparing the original network with the individual layer added network
individual_layer_plot_physical <- generate_combined_stage_plot(original_cov_data_physical, individual_cov_data_physical, "Individual Layer Links Added - Physical")
print(individual_layer_plot_physical)

# Analyze and plot network after adding two layer combination links
two_layer_cov_data_physical <- process_network(extended_edge_lists_two_layers, "Two Layer - Directed - Physical", num_layers, "physical", TRUE)
two_layer_combined_plot_physical <- generate_combined_stage_plot(original_cov_data_physical, two_layer_cov_data_physical, "Two Layer Combination Links Added")
print(two_layer_combined_plot_physical)

#Analyzing and plotting network after adding three layer combination links
three_layer_cov_data_physical <- process_network(extended_edge_lists_three_layers, "Three Layer Combination Links Added", num_layers, "physical", TRUE)
three_layer_combined_plot_physical <- generate_combined_stage_plot(original_cov_data_physical, three_layer_cov_data_physical, "Three Layer Combination Links Added")
print(three_layer_combined_plot_physical)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Relaxed Physical ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

# Analyze and visualize the original network (Directed, Relaxed Physical)
original_cov_data_relaxed_physical <- process_network(extended_edge_lists_original, "Original - Directed - Relaxed Physical", num_layers, "relaxed-physical", FALSE)

# Process the network with individual layer edges added (Directed, Relaxed Physical)
individual_cov_data_relaxed_physical <- process_network(extended_edge_lists_individual, "Individual Layer - Directed - Relaxed Physical", num_layers, "relaxed-physical", TRUE)

# Generate and print the plot comparing the original network with the individual layer added network
individual_layer_plot_relaxed_physical <- generate_combined_stage_plot(original_cov_data_relaxed_physical, individual_cov_data_relaxed_physical, "Individual Layer Links Added - Relaxed Physical")
print(individual_layer_plot_relaxed_physical)

# Analyze and plot network after adding two layer combination links
two_layer_cov_data_relaxed_physical <- process_network(extended_edge_lists_two_layers, "Two Layer - Directed - Relaxed Physical", num_layers, "relaxed-physical", TRUE)
two_layer_combined_plot_relaxed_physical <- generate_combined_stage_plot(original_cov_data_relaxed_physical, two_layer_cov_data_relaxed_physical, "Two Layer Combination Links Added")
print(two_layer_combined_plot_relaxed_physical)

#Analyzing and plotting network after adding three layer combination links
three_layer_cov_data_relaxed_physical <- process_network(extended_edge_lists_three_layers, "Three Layer Combination Links Added", num_layers, "relaxed-physical", TRUE)
three_layer_combined_plot_relaxed_physical <- generate_combined_stage_plot(original_cov_data_relaxed_physical, three_layer_cov_data_relaxed_physical, "Three Layer Combination Links Added")
print(three_layer_combined_plot_relaxed_physical)


# --------------------------------------------------------------------------------------------------------------
                          # Undirected Network Analysis Using Other RW Strategies
# --------------------------------------------------------------------------------------------------------------


############# Step 05 Final Combined Plot#############

# Final combined plot of all stages
all_stages_combined_data <- rbind(
  transform(original_cov_data_pagerank, Stage = "Original"),
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
final_combined_plot <- generate_final_combined_plot(original_cov_data_pagerank, individual_cov_data, two_layer_cov_data, three_layer_cov_data)
print(final_combined_plot)
