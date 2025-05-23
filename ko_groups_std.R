library(dplyr)
library(tidyr)
library(ggplot2)
library(dbscan)
library(umap)
library(colorspace)
library(uwot)


# Load data
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

# Standardize the data
binary_matrix_std <- scale(binary_matrix_filtered, center = TRUE, scale = FALSE) 

set.seed(42)

umap_res <- umap(binary_matrix_std, n_neighbors = 15, min_dist = 0.3, metric = "hamming")
umap_df <- data.frame(Plasmid = rownames(binary_matrix_filtered),
                      UMAP1 = umap_res$layout[,1],
                      UMAP2 = umap_res$layout[,2])

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "UMAP of Plasmid KO Profiles")


# DBSCAN Clustering on UMAP
dbscan_res <- dbscan(umap_df[,2:3], eps = 1, minPts = 8)
umap_df$Cluster <- as.factor(dbscan_res$cluster)

# Compute centroids for each cluster
cluster_centroids <- umap_df |> 
  group_by(Cluster) |> 
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

# Plot UMAP with DBSCAN clusters
umap_cluster_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(Cluster))) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text(data = cluster_centroids, aes(label = Cluster), 
            color = "black", fontface = "bold", size = 4) +
  labs(title = "DBSCAN Clustering on Standardized UMAP", color = "Cluster") +
  theme_minimal()

ggsave("C:/Users/hayat/Downloads/R_files/graphs/UMAP_DBSCAN_std.png", plot = umap_cluster_plot, width = 7, height = 5, dpi = 300, bg = "white")

# Save clustering results
write.csv(umap_df, "C:/Users/hayat/Downloads/R_files/data/UMAP_DBSCAN_clusters_std.csv", row.names = FALSE)

# Generate a list of plasmids per cluster
plasmid_list <- split(umap_df$Plasmid, umap_df$Cluster)

# Write to a text file
sink("C:/Users/hayat/Downloads/R_files/data/plasmids_per_cluster_std.txt")
for (cluster in names(plasmid_list)) {
  cat(paste0("Cluster ", cluster, ":\n"))
  cat(paste(plasmid_list[[cluster]], collapse = "\n"), "\n\n")
}
sink()


any(is.na(binary_matrix_std))  # Should return FALSE
any(is.infinite(binary_matrix_std))  # Should return FALSE

summary(umap_res$layout)
plot(umap_res$layout[,1], umap_res$layout[,2])  # Quick check before ggplot

