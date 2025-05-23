# Load necessary libraries
library(ggplot2)

# Load your data (adjust the path and specify tab separator for TSV)
df <- read.table("data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv", sep = "\t", header = TRUE, row.names = 1)

# Remove NAs
df <- na.omit(df)

# Assuming your data is stored in a variable 'data'
# Remove any non-numeric columns if needed
data_numeric <- df[, sapply(df, is.numeric)]

# Identify columns with zero variance
zero_variance_cols <- sapply(data_numeric, function(x) var(x) == 0)

# Remove columns with zero variance
data_numeric <- data_numeric[, !zero_variance_cols]

# Add pseudocount for better handling of zero values
data_numeric[data_numeric == 0] <- 1e-6

# Standardize the data (optional but recommended)
data_scaled <- scale(data_numeric)

# Remove sample names (row names) from the data
data_scaled <- as.data.frame(data_scaled)

# Perform PCA
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# Check the summary of the PCA
summary(pca_result)

# Calculating the variance for each PC
variance <- pca_result$sdev^2

# Calculating the proportion of variance explained by each component
explained_variance <- variance / sum(variance)

# Summary table
summary_table <- data.frame(PC = 1:length(explained_variance),
                            Variance = variance,
                            Explained_Variance = explained_variance,
                            Cumulative_Variance = cumsum(explained_variance))

print(summary_table)


# Scree plot to visualize the variance explained by each PC
screeplot(pca_result, main = "Scree Plot", col = "blue", type = "lines")

png("graphs/PCA_biplot.png", width = 800, height = 600)  

# Biplot (no sample names shown)
biplot(pca_result, main = "PCA Biplot", cex = 0.1, col = c("blue", "red"))


dev.off()

# Add scatter plot of first two principal components (PC1 and PC2)
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Sample = rownames(df))
png("graphs/PCA_scatterplot.png", width = 800, height = 600)
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(color = "blue") +
  geom_text(aes(label = ifelse(abs(PC1) > 0.5 | abs(PC2) > 0.5, as.character(Sample), "")), 
            hjust = 1.1, vjust = 1.1, size = 3) + 
  labs(title = "PCA Scatter Plot", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

dev.off()

# Calculate proportions
rpkm_proportions <- sweep(data_numeric, 2, colSums(data_numeric), FUN = "/")

# Perform PCA on this proportion matrix without scaling
pca_result <- prcomp(rpkm_proportions, scale = FALSE)

png("graphs/PCA_proportions.png", width = 800, height = 600)
# To visualize the PCA
plot(pca_result$x[, 1:2], main = "PCA of RPKM Proportions", xlab = "PC1", ylab = "PC2")

dev.off()

col_sums <- colSums(data_numeric)
summary(col_sums)

# Plot PC3 vs PC4
png("graphs/PCA_proportions_3_4.png", width = 800, height = 600)
plot(pca_result$x[, 3], pca_result$x[, 4], xlab = "PC3", ylab = "PC4")
dev.off()


# Calculate min, max, and mean (or median) for each column
summary_stats <- data.frame(
  Min = apply(rpkm_proportions, 2, min),
  Max = apply(rpkm_proportions, 2, max),
  Mean = apply(rpkm_proportions, 2, mean),      # Use median() if preferred
  Median = apply(rpkm_proportions, 2, median)
)

# Print the summary statistics
print(summary_stats)

# Log transformation (adding a small constant to avoid log(0))
data_log <- log(data_numeric)  # Adjust constant if needed

# Perform PCA on the log-transformed data
pca <- prcomp(data_log, center = TRUE, scale. = TRUE)

# Summary of PCA to check the variance explained by each principal component
summary(pca)

# Create a data frame with the PCA results
pca_df <- data.frame(pca$x)  # Principal component scores
png("graphs/PCA_log.png", width = 800, height = 600)
# Plot the first two principal components
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "blue", alpha = 0.7) +
  theme_minimal() +
  ggtitle("PCA Plot: PC1 vs. PC2") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2")
dev.off()
# Variance explained by each component (scree plot)

png("graphs/scree_log.png", width = 800, height = 600)
ggplot(data.frame(PC = 1:length(pca$sdev), Variance = pca$sdev^2 / sum(pca$sdev^2) * 100), 
       aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  theme_minimal() +
  ggtitle("Scree Plot: Variance Explained by Each PC") +
  xlab("Principal Component") +
  ylab("Variance (%)")

dev.off()


