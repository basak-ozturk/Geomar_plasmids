# Set CRAN mirror at the beginning of the script
#options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load necessary packages
#install.packages("pheatmap")  # Run only once if not installed
#install.packages("gplots")

# Load required libraries
library(Matrix)
library(pheatmap)

# Load your data (adjust the path and specify tab separator for TSV)
df <- read.table("data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv", sep = "\t", header = TRUE, row.names = 1)

# Check the dimensions before omitting NAs
cat("Original data dimensions:", dim(df), "\n")

# Remove NAs
df <- na.omit(df)


# Log-transform the data to reduce zeros' effect
df_log <- log1p(df)

df_scaled <- scale(df_log)



# # Convert to sparse matrix
# df_sparse <- as(df_scaled, "dgCMatrix")
# 
# # Convert back to matrix for pheatmap (efficient memory use)
# df_sparse_matrix <- as.matrix(df_sparse)
# 
# # Remove row and column names to save memory
# rownames(df_sparse_matrix) <- NULL
# colnames(df_sparse_matrix) <- NULL
# 
# # Generate heatmap without labels
# pheatmap(df_sparse_matrix, show_rownames = FALSE, show_colnames = FALSE)

write.csv(df_scaled[1:500, 1:50], "data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput_sample.tsv")


