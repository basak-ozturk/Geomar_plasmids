# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster)
library(vegan)
library(dbscan)
library(Rtsne)
library(umap)
library(colorspace)

df <- read.table("C:/Users/hayat/Downloads/R_files/data/KEGG_ko_per_plasmid.txt", 
                 header=FALSE, sep="\t", col.names=c("Plasmid", "KO"))

# Split KO into a list
df <- df |> mutate(KO = strsplit(KO, ","))

# Create a binary presence/absence matrix
binary_matrix <- df |> 
  unnest(KO) |> 
  mutate(Value = 1) |> 
  pivot_wider(names_from = KO, values_from = Value, values_fill = 0)

# Convert to matrix format
rownames(binary_matrix) <- binary_matrix$Plasmid
binary_matrix <- binary_matrix[, -1]

# Filtering KO that appear in at least 10 plasmids
filtered_columns <- colSums(binary_matrix) >= 10
binary_matrix_filtered <- binary_matrix[, filtered_columns]

# Perform PCA on filtered data
pca_res_filtered <- prcomp(binary_matrix_filtered, scale. = FALSE)
explained_var <- summary(pca_res_filtered)$importance[2,] * 100


# K-Means clustering
set.seed(42)
kmeans_res <- kmeans(binary_matrix_filtered, centers = 3, nstart = 25)
binary_matrix_filtered$Cluster <- as.factor(kmeans_res$cluster)

# Create a dataframe for visualization
pca_df_filtered <- data.frame(
  Plasmid = rownames(binary_matrix_filtered),
  PC1 = pca_res_filtered$x[,1],
  PC2 = pca_res_filtered$x[,2],
  Cluster = as.factor(kmeans_res$cluster)  # Use factors for discrete clusters
)

# PCA plot with k-means clusters
pca_plot <- ggplot(pca_df_filtered, aes(x = PC1, y = PC2, color = Cluster, label = Plasmid)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text(hjust = 0.5, vjust = -0.5, size = 3) +
  labs(title = "PCA of Plasmid KO Clusters (Filtered & K-Means)",
       x = paste0("PC1 (", round(explained_var[1], 2), "%)"),
       y = paste0("PC2 (", round(explained_var[2], 2), "%)"),
       color = "Cluster") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")  # Adjust colors for better visibility

# Save the PCA plot
ggsave("C:/Users/hayat/Downloads/R_files/graphs/PCA_KMeans_plot_ko.png", plot = pca_plot, width = 5, height = 4, dpi = 300, bg = "white")


# UMAP dimensionality reduction
set.seed(42)
binary_matrix_filtered_numeric <- as.matrix(binary_matrix_filtered)
mode(binary_matrix_filtered_numeric) <- "numeric"
umap_res <- umap(binary_matrix_filtered_numeric, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

# Convert UMAP results into a dataframe
umap_df <- data.frame(Plasmid = rownames(binary_matrix_filtered),
                      UMAP1 = umap_res$layout[,1],
                      UMAP2 = umap_res$layout[,2])

# Plot UMAP results
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 1, alpha = 0.7) +
  labs(title = "UMAP Visualization of Plasmid KO") +
  theme_minimal()

# DBSCAN Clustering on UMAP
dbscan_res <- dbscan(umap_df[,2:3], eps = 1, minPts = 8)
umap_df$Cluster <- as.factor(dbscan_res$cluster)

cluster_colors <- qualitative_hcl(33, palette = "Dark 3")
cluster_centroids <- umap_df |> # Compute centroids for each cluster
  group_by(Cluster) |>
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

# Plot UMAP with cluster labels
umap_plot<-ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(Cluster))) +
  geom_point(size = 1, alpha = 0.7) +  # Use alpha to make overlapping points more visible
  geom_text(data = cluster_centroids, aes(label = Cluster), 
            color = "black", fontface = "bold", size = 4) +  # Cluster number labels
  labs(title = "DBSCAN Clustering of UMAP Results", color = "Cluster") +
  theme_minimal() +
  scale_color_manual(values = cluster_colors) +
  theme(legend.position = "right")


# Save with larger size
ggsave("C:/Users/hayat/Downloads/R_files/graphs/UMAP_DBSCAN_plot_ko.png", plot = umap_plot, width = 7, height = 5, dpi = 300, bg = "white")

densityplot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.5) +  # Scatter plot with transparency
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) +
  scale_fill_viridis_c() +  # Use a color gradient
  theme_minimal() +
  labs(title = "UMAP Density Visualization")

# Save with larger size
ggsave("C:/Users/hayat/Downloads/R_files/graphs/UMAP_DBSCAN_density_plot_ko.png", plot = densityplot, width = 7, height = 5, dpi = 300, bg = "white")


# Step 3: Save clustering results
write.csv(umap_df, "C:/Users/hayat/Downloads/R_files/data/UMAP_DBSCAN_clusters_ko.csv", row.names = FALSE)

# Step 4: Now read the file (only if needed)
umap_results <- read.csv("C:/Users/hayat/Downloads/R_files/data/UMAP_DBSCAN_clusters_ko.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Step 5: Ensure correct Plasmid assignment
umap_results$Plasmid <- df$Plasmid

# Step 6: Reorder columns
umap_results <- umap_results[, c("Plasmid", "UMAP1", "UMAP2", "Cluster")]

# Step 7: Save again
write.csv(umap_results, "C:/Users/hayat/Downloads/R_files/data/UMAP_DBSCAN_clusters_ko.csv", row.names = FALSE)



umap_results <- read.csv("C:/Users/hayat/Downloads/R_files/data/UMAP_DBSCAN_clusters_ko.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Ensure the row numbers in UMAP match the order in the original plasmid list
umap_results$Plasmid <- df$Plasmid

# Reorder columns to make it clearer
umap_results <- umap_results[, c("Plasmid", "UMAP1", "UMAP2", "Cluster")]

# Save clustering results
write.csv(umap_results, "C:/Users/hayat/Downloads/R_files/data/UMAP_DBSCAN_clusters_ko.csv", row.names = FALSE)


#Get a list of plasmids for each cluster
plasmid_list <-split(umap_results$Plasmid, umap_results$Cluster)

# Write to a text file
sink("C:/Users/hayat/Downloads/R_files/data/plasmids_per_cluster_ko.txt")
for (cluster in names(plasmid_list)) {
  cat(paste0("Cluster ", cluster, ":\n"))
  cat(paste(plasmid_list[[cluster]], collapse = "\n"), "\n\n")
}
sink()  # Close the file connection




