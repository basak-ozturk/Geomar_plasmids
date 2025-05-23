import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

# Load data
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_pathway_per_plasmid_with_names.txt", sep="\t", header=None, names=["Plasmid", "Pathway"])

# Convert Pathway to list
df["Pathway"] = df["Pathway"].apply(lambda x: x.split(","))

# Create binary presence/absence matrix
binary_matrix = df.explode("Pathway").assign(Value=1).pivot(index="Plasmid", columns="Pathway", values="Value").fillna(0)

# Filter K-numbers that appear in at least 100 plasmids
filtered_columns = binary_matrix.columns[binary_matrix.sum(axis=0) >= 10]
binary_matrix_filtered = binary_matrix[filtered_columns]

# Compute distance matrix (Jaccard distance)
distance_matrix = pdist(binary_matrix_filtered, metric="euclidean")

# Perform hierarchical clustering (change linkage method if needed)
linkage_matrix = linkage(distance_matrix, method="ward")  # Use 'average' for Jaccard

# Define a cut-off threshold for clustering (adjust as needed)
cluster_cutoff = 0.5  # Adjust based on the dendrogram structure

# Assign cluster labels
cluster_labels = fcluster(linkage_matrix, t=cluster_cutoff, criterion="distance")

# Add cluster labels to dataframe
binary_matrix_filtered["Cluster"] = cluster_labels

# Save cluster assignments
binary_matrix_filtered.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_clusters_hierarchical.csv")

# Set Seaborn style
sns.set_style("whitegrid")

# Plot dendrogram and save it
plt.figure(figsize=(12, 6))
dendrogram(linkage_matrix, no_labels=True)  # Show labels if needed
plt.axhline(y=cluster_cutoff, color='r', linestyle='--')  # Show cut-off threshold
plt.title("Hierarchical Clustering Dendrogram Based on Encoded Pathways")
plt.xlabel("Plasmids")
plt.ylabel("Distance")
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/pathway_dendrogram_10.png", dpi=300, bbox_inches="tight")  # Save as PNG
plt.show()

clustermap = sns.clustermap(
    binary_matrix_filtered.drop(columns=["Cluster"]),
    row_cluster=True,
    col_cluster=False,
    cmap="coolwarm",
    figsize=(12, 8),
    xticklabels=False,
    yticklabels=False,
    dendrogram_ratio=(0.1, 0.2),  # Reduce dendrogram size
    cbar_pos=(1, 0.3, 0.03, 0.5)  # Move colorbar down
)

# Rotate x-axis labels
plt.setp(clustermap.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

# Save figure
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/pathway_clustermap.png", dpi=300, bbox_inches="tight")  
plt.show()


# Save plasmid names grouped by cluster (FIXED)
cluster_groups = binary_matrix_filtered.reset_index().groupby("Cluster")["Plasmid"].apply(lambda x: ", ".join(x))

# Save to a CSV file
cluster_groups.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_10.csv", header=True)

print("Plasmid names per cluster saved to 'plasmid_cluster_mapping_10.csv'")

# Load the cluster mapping file
cluster_mapping = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_10.csv", 
    index_col=0, 
    names=["Cluster", "Plasmids"], 
    skiprows=1  # Adjust if needed
)

# Handle missing values in "Plasmids" column
cluster_mapping["Plasmids"] = cluster_mapping["Plasmids"].fillna("")

# Count the number of plasmids per cluster
cluster_mapping["Plasmid_Count"] = cluster_mapping["Plasmids"].apply(lambda x: len(str(x).split(", ")) if x else 0)

# Get total number of plasmids
total_plasmids = cluster_mapping["Plasmid_Count"].sum()
print(f"Total number of plasmids: {total_plasmids}")

# Calculate percentage representation per cluster
cluster_mapping["Percentage"] = (cluster_mapping["Plasmid_Count"] / total_plasmids) * 100

# Filter clusters with more than 50 plasmids
filtered_clusters = cluster_mapping[cluster_mapping["Plasmid_Count"] > 50]

# Save filtered clusters
filtered_clusters.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_filtered_10.csv", 
    index=True,  
    header=True,  
    columns=["Plasmid_Count", "Percentage"]
)

print("\nFiltered clusters (more than 10 plasmids) saved to 'plasmid_cluster_filtered.csv'.")

# Sort clusters by cluster number for a clearer x-axis
filtered_clusters_sorted = filtered_clusters.sort_index()

# Plot bar chart (Cluster Number on x-axis, Plasmid Count on y-axis)
plt.figure(figsize=(12, 6))
sns.barplot(
    x=filtered_clusters_sorted.index.astype(str),  # Convert to string for better labels
    y=filtered_clusters_sorted["Plasmid_Count"],
    color="blue",
    edgecolor="black"
)

# Save histogram
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_cluster_barplot_10.png", dpi=300, bbox_inches="tight")
plt.show()

