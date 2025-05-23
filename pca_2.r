
library(ggplot2)
library(caret)


df <- read.table("data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv", sep = "\t", header = TRUE, row.names = 1)

# Remove NAs
df <- na.omit(df)

# Remove any non-numeric columns if needed
data_numeric <- df[, sapply(df, is.numeric)]

# Identify columns with zero variance
zero_variance_cols <- sapply(data_numeric, function(x) var(x) == 0)

# Remove columns with zero variance
data_numeric <- data_numeric[, !zero_variance_cols]

# Add pseudocount for better handling of zero values
data_numeric[data_numeric == 0] <- 1e-6

# **Transpose first, then scale (correct order)**
data_scaled <- scale(t(data_numeric))

# Perform PCA on scaled data
pca_result <- prcomp(data_scaled, center = TRUE, scale. = FALSE)  # Don't scale again inside `prcomp`

# Convert PCA results into a data frame
pca_df <- as.data.frame(pca_result$x)
pca_df$Total_Abundance <- total_abundance  # Add abundance as a column

png("graphs/PCA_RPKM_abundance.png", width = 800, height = 600)  

# PCA Plot with coloring by total RPKM abundance
ggplot(pca_df, aes(x = PC1, y = PC2, color = Total_Abundance)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "PCA of Metagenome RPKM Data (Scaled)",
       x = "PC1",
       y = "PC2",
       color = "Total RPKM Abundance") +
  theme_minimal()

dev.off()

# Check the loadings (contribution of each feature to each principal component)
head(pca_result$rotation)


# Min-Max Scaling
preproc <- preProcess(data_numeric, method = "range")
min_max_scaled_data <- predict(preproc, data_numeric)

# Log Transformation (if appropriate)
log_scaled_data <- log(data_numeric)

# Compare PCA on each scaled dataset
pca_standardized <- prcomp(data_scaled)  # Using  already standardized data
pca_min_max <- prcomp(min_max_scaled_data)
pca_log <- prcomp(log_scaled_data)

# Plot PCA results for comparison
par(mfrow = c(1, 3))  # Arrange multiple plots in one row
plot(pca_standardized$x[, 1:2], main = "Standardized PCA")
plot(pca_min_max$x[, 1:2], main = "Min-Max Scaled PCA")
plot(pca_log$x[, 1:2], main = "Log Transformed PCA")


set.seed(123)  # For reproducibility
kmeans_result <- kmeans(pca_result$x[, 1:2], centers = 4)  # Adjust 'centers' as needed

# Combine PCA results and clustering labels into a data frame
pca_result_with_clusters <- data.frame(pca_result$x[, 1:2], cluster = factor(kmeans_result$cluster))

png("graphs/PCA_K_means_4.png", width = 800, height = 600)  

# Plot the PCA results with clustering
plot(pca_result_with_clusters[, 1:2], col = pca_result_with_clusters$cluster, pch = 19,
     main = "PCA with K-means Clustering", xlab = "PC1", ylab = "PC2", cex=1.5)

dev.off()


# Calculate the explained variance for each PC
explained_standardized <- (pca_standardized$sdev^2) / sum(pca_standardized$sdev^2) * 100
explained_min_max <- (pca_min_max$sdev^2) / sum(pca_min_max$sdev^2) * 100
explained_log <- (pca_log$sdev^2) / sum(pca_log$sdev^2) * 100

# Print variance explained for PC1 and PC2
cat("Variance Explained by PC1 and PC2 (Standardized Data):\n")
cat("PC1: ", explained_standardized[1], "%\n")
cat("PC2: ", explained_standardized[2], "%\n\n")

cat("Variance Explained by PC1 and PC2 (Min-Max Scaled Data):\n")
cat("PC1: ", explained_min_max[1], "%\n")
cat("PC2: ", explained_min_max[2], "%\n\n")

cat("Variance Explained by PC1 and PC2 (Log Transformed Data):\n")
cat("PC1: ", explained_log[1], "%\n")
cat("PC2: ", explained_log[2], "%\n\n")

# Plot PCA results for comparison
par(mfrow = c(1, 3))  # Arrange multiple plots in one row
plot(pca_standardized$x[, 1:2], main = paste("Standardized PCA\nPC1: ", round(explained_standardized[1], 2), "%, PC2: ", round(explained_standardized[2], 2), "%"))
plot(pca_min_max$x[, 1:2], main = paste("Min-Max Scaled PCA\nPC1: ", round(explained_min_max[1], 2), "%, PC2: ", round(explained_min_max[2], 2), "%"))
plot(pca_log$x[, 1:2], main = paste("Log Transformed PCA\nPC1: ", round(explained_log[1], 2), "%, PC2: ", round(explained_log[2], 2), "%"))

# Create scree plot for Standardized Data
png("graphs/scree_standardized.png", width = 800, height = 600)
plot(explained_standardized, type = "b", main = "Scree Plot (Standardized Data)", 
     xlab = "Principal Components", ylab = "Variance Explained (%)", pch = 19, col = "blue", cex = 1.5)
dev.off()

# Create scree plot for Min-Max Scaled Data
png("graphs/scree_min_max.png", width = 800, height = 600)
plot(explained_min_max, type = "b", main = "Scree Plot (Min-Max Scaled Data)", 
     xlab = "Principal Components", ylab = "Variance Explained (%)", pch = 19, col = "green", cex = 1.5)
dev.off()

# Create scree plot for Log Transformed Data
png("graphs/scree_log.png", width = 800, height = 600)
plot(explained_log, type = "b", main = "Scree Plot (Log Transformed Data)", 
     xlab = "Principal Components", ylab = "Variance Explained (%)", pch = 19, col = "red", cex = 1.5)
dev.off()



set.seed(123)  # For reproducibility
kmeans_log <- kmeans(pca_log$x[, 1:2], centers = 3)  # Adjust 'centers' as needed

# Combine PCA results and clustering labels into a data frame
pca_log_with_clusters <- data.frame(pca_log$x[, 1:2], cluster = factor(kmeans_log$cluster))  # Corrected here

# Save the plot to a file
png("graphs/PCA_log_K_means_3.png", width = 800, height = 600)  

# Plot the PCA results with clustering
plot(pca_log_with_clusters[, 1:2], col = pca_log_with_clusters$cluster, pch = 19,
     main = "PCA with K-means Clustering", xlab = "PC1", ylab = "PC2", cex = 1.5)

# Close the graphics device
dev.off()

