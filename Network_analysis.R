library(corrplot)
library(igraph)
library(RColorBrewer)

Over_1000 <- read.csv("PHQ9_cluster_networks.csv")


All_net <- readRDS("All_GTEx_DLPFC_networks.RData")
df_names <- names(All_net)
new_df_names <- sub("_+$", "", df_names)
names(All_net) <- new_df_names

process_dataframe <- function(df, df_name) {
  df <- df[, c("Gene", "R"), drop = FALSE]  # Extract "Gene" and "R" columns
  colnames(df)[colnames(df) == "R"] <- df_name  # Rename the "R" column to dataframe name
  return(df)
}

# Applying the function to all dataframes in the list
list_of_dataframes <- lapply(names(All_net), function(name){
  process_dataframe(All_net[[name]], name)
})

data <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), list_of_dataframes)

data <- column_to_rownames(data, var = "Gene")

data <- data[,colnames(data) %in% Over_1000$ensembl_gene_id]


# Make correlation matrix
correlation_matrix <- round(cor(data),2)


# Convert correlation matrix to adjacency matrix
adjacency_matrix_pos_neg <- as.matrix(correlation_matrix)
adjacency_matrix_pos_neg <- ifelse(adjacency_matrix_pos_neg < 0.8 & adjacency_matrix_pos_neg > -0.8, 0, adjacency_matrix_pos_neg)
adjacency_matrix <- abs(adjacency_matrix_pos_neg)

# Create graph object
network <- graph.adjacency(adjacency_matrix, mode = "upper", weighted = TRUE, diag = FALSE)
non_network <- graph.adjacency(adjacency_matrix_pos_neg, mode = "upper", weighted = TRUE, diag = FALSE)

# Define the range for edge thickness
min_thickness <- 0.5  # Minimum edge thickness
max_thickness <- 5    # Maximum edge thickness

# Normalize edge weights to the range [min_thickness, max_thickness]
normalized_weights <- E(network)$weight
normalized_weights <- (normalized_weights - min(normalized_weights)) / 
  (max(normalized_weights) - min(normalized_weights))
edge_thickness <- min_thickness + normalized_weights * (max_thickness - min_thickness) 


# Define edge colors based on positive/negative correlation
edge_colors <- ifelse(E(non_network)$weight < 0, "skyblue", "#FF7F7F")

node_size <- data %>% 
  summarise(across(everything(), ~ sum(. > 0.5 | . < -0.5, na.rm = TRUE)))
min_thickness <- 6  # Minimum edge thickness
max_thickness <- 10 
node_size <- (node_size - min(node_size)) / 
  (max(node_size) - min(node_size))
node_size <- min_thickness + node_size * (max_thickness - min_thickness)
node_size <- t(node_size)

# Reorder the rows in 'Over_1000' based on the order of 'data_colnames'
data_colnames <- colnames(data)
Over_1000 <- Over_1000[order(match(Over_1000$ensembl_gene_id, data_colnames)), ]

# Extract node names and cluster information from the Networks dataframe
color_palette <- c("palegreen3","plum3","gold2")
Over_1000$Cluster <- factor(Over_1000$Cluster, levels = unique(Over_1000$Cluster))
node_colors <- color_palette[Over_1000$Cluster]

# Label FADS1
labels <- rep(NA, vcount(network))
labels[13] <- "FADS1 network"


betweenness <- betweenness(network) %>% 
  as.data.frame()%>%
  dplyr::rename(Betweenness='.')%>%
  rownames_to_column() 


