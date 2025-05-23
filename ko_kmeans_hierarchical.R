library(tidyverse)
library(cluster)
library(ggdendro)
library(ggplot2)
library(pheatmap)
library(cluster)
library(dendextend)
library(dplyr)
# Load Data
df <- read_delim("C:/Users/hayat/Downloads/R_files/data/KEGG_pathway_per_plasmid_with_names.txt", 
                 delim = "\t", col_names = c("Plasmid", "Pathway"))

# Convert Pathway column to list
df <- df |> 
  mutate(Pathway = str_split(Pathway, ",")) 

# Create a binary presence/absence matrix
binary_matrix <- df |> 
  unnest(Pathway) |> 
  mutate(Value = 1) |> 
  pivot_wider(names_from = Pathway, values_from = Value, values_fill = 0)

# Preserve Plasmid as rownames
binary_matrix <- binary_matrix |> 
  column_to_rownames(var = "Plasmid")

# Filter pathways that appear in at least 100 plasmids
filtered_matrix <- binary_matrix[, colSums(binary_matrix) >= 100]

# Ensure filtered_matrix is a data frame and reset rownames as a column
filtered_matrix <- as.data.frame(filtered_matrix) |> 
  rownames_to_column(var = "Plasmid")

# Ensure only numeric columns are used for clustering
filtered_numeric <- filtered_matrix |> select(-Plasmid)

# Check for missing values before K-means
if (any(is.na(filtered_numeric))) {
  stop("Error: NA values found in filtered_numeric matrix.")
}



silhouette_scores <- vector()
for (k in 2:20) {  # K must be at least 2
  kmeans_result <- kmeans(filtered_numeric, centers = k, nstart = 25, iter.max = 20)
  sil <- silhouette(kmeans_result$cluster, dist(filtered_numeric))
  silhouette_scores[k] <- mean(sil[, 3])
}

best_k <- which.max(silhouette_scores)  # Finds the k with the highest score
print(paste("Best K:", best_k))


# Plot silhouette scores
plot(2:20, silhouette_scores[-1], type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (K)", ylab = "Average Silhouette Score",
     main = "Silhouette Method for Optimal K")


# Perform K-means clustering
set.seed(123)  # For reproducibility
num_clusters <- 4  # Adjust as needed
kmeans_result <- kmeans(filtered_numeric, centers = num_clusters)

# Add K-means cluster labels to filtered_matrix
filtered_matrix$Kmeans_Cluster <- as.factor(kmeans_result$cluster)

# Debug: Check if cluster labels exist
table(filtered_matrix$Kmeans_Cluster)

# Standardize Plasmid names for matching
dend_labels <- tolower(trimws(filtered_matrix$Plasmid))
filtered_matrix$Plasmid <- dend_labels
rownames(filtered_numeric) <- dend_labels  # Ensure rownames match

# Perform hierarchical clustering
hc <- hclust(dist(filtered_numeric, method = "manhattan"))  # Try "manhattan" instead of default "euclidean"


# Convert hclust object to dendrogram
dend <- as.dendrogram(hc)
dend_data <- ggdendro::dendro_data(dend)

# Standardize dendrogram labels for matching
dend_data$labels$label <- tolower(trimws(dend_data$labels$label))

# Assign correct cluster labels
dend_data$labels <- dend_data$labels |> 
  mutate(Kmeans_Cluster = filtered_matrix$Kmeans_Cluster[match(label, filtered_matrix$Plasmid)])

# Check if any clusters are still NA
if (any(is.na(dend_data$labels$Kmeans_Cluster))) {
  stop("Error: Mismatch between dendrogram labels and filtered_matrix plasmid names.")
}

filtered_matrix$Kmeans_Cluster <- as.factor(filtered_matrix$Kmeans_Cluster)
# Ensure cluster colors are assigned to existing levels only
cluster_levels <- levels(filtered_matrix$Kmeans_Cluster)

cluster_colors <- setNames(rainbow(length(levels(filtered_matrix$Kmeans_Cluster))), 
                           levels(filtered_matrix$Kmeans_Cluster))

ggplot() +
  geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  
  # Add solid color blocks instead of individual colored text
  geom_tile(data = dend_data$labels, aes(x = x, y = -0.5, fill = as.factor(Kmeans_Cluster)), height = 0.5) +
  
  # Keep plasmid labels black for readability
  #geom_text(data = dend_data$labels, aes(x = x, y = y, label = label), angle = 90, hjust = 1, size = 3, color = "black") +
  
  scale_fill_manual(values = cluster_colors, 
                    guide = guide_legend(override.aes = list(shape = 21, size = 5))) +
  
  theme_minimal() +
  labs(title = "Hierarchical Clustering of Plasmids based on ko numbers with K-means Clusters",
       x = "Plasmid", y = "Distance", fill = "K-means Cluster") +
  
  # Remove x-axis tick labels
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



# Save dendrogram plot
ggsave("C:/Users/hayat/Downloads/R_files/graphs/pathway_dendrogram_kmeans_14.png", dpi = 300, width = 12, height = 6, bg="white")



# Prepare data for pheatmap
heatmap_data <- as.matrix(filtered_matrix |> select(-Plasmid, -Kmeans_Cluster))
rownames(heatmap_data) <- filtered_matrix$Plasmid

# Prepare annotation for pheatmap
annotation_df <- data.frame(Kmeans_Cluster = filtered_matrix$Kmeans_Cluster)
rownames(annotation_df) <- filtered_matrix$Plasmid

# Define annotation colors
ann_colors <- list(Kmeans_Cluster = cluster_colors)

# Save heatmap
png("C:/Users/hayat/Downloads/R_files/graphs/pathway_clustermap_kmeans.png", width = 1200, height = 800, res = 300)
pheatmap(heatmap_data, 
         cluster_rows = hc, 
         cluster_cols = FALSE, 
         annotation_row = annotation_df,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# Save Plasmid Clusters
write_csv(filtered_matrix, "C:/Users/hayat/Downloads/R_files/data/plasmid_clusters_kmeans.csv")
