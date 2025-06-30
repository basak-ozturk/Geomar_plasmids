import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
import numpy as np
#import matplotlib.colors as mcolors

# Load data
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_pathway_per_plasmid_with_names.txt", sep="\t", header=None, names=["Plasmid", "Pathway"])

# Convert Pathway to list
df["Pathway"] = df["Pathway"].apply(lambda x: x.split(","))

# Create binary presence/absence matrix
binary_matrix = df.explode("Pathway").assign(Value=1).pivot(index="Plasmid", columns="Pathway", values="Value").fillna(0)

# Filter K-numbers that appear in at least 100 plasmids
filtered_columns = binary_matrix.columns[binary_matrix.sum(axis=0) >= 100]
binary_matrix_filtered = binary_matrix[filtered_columns]

# Compute distance matrix (Jaccard distance can be used, too)
distance_matrix = pdist(binary_matrix_filtered, metric="euclidean")

# Perform hierarchical clustering (change linkage method if needed)
linkage_matrix = linkage(distance_matrix, method="ward")  # Use 'average' for Jaccard

# Define a cut-off threshold for clustering (adjust as needed)
cluster_cutoff = 26  # Adjust based on the dendrogram structure

# Assign cluster labels
cluster_labels = fcluster(linkage_matrix, t=cluster_cutoff, criterion="distance")

# Add cluster labels to dataframe
binary_matrix_filtered["Cluster"] = cluster_labels

# Save cluster assignments
binary_matrix_filtered.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_clusters_hierarchical_100.csv")

# Set Seaborn style
sns.set_style("whitegrid")

# Plot dendrogram and save it
plt.figure(figsize=(12, 6))
#dendrogram(linkage_matrix, no_labels=True)  # Show labels if needed
dendrogram(linkage_matrix, no_labels=True, color_threshold=cluster_cutoff)
plt.axhline(y=cluster_cutoff, color='r', linestyle='--')  # Show cut-off threshold
plt.title("Hierarchical Clustering Dendrogram Based on Encoded Pathways")
plt.xlabel("Plasmids")
plt.ylabel("Distance")
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/pathway_dendrogram_100.png", dpi=300, bbox_inches="tight")  # Save as PNG
plt.show()

clustermap = sns.clustermap(
    binary_matrix_filtered.drop(columns=["Cluster"]),
    row_cluster=True,   # Clustering rows
    col_cluster=True,   # Clustering columns as well
    cmap="Greys",
    figsize=(12, 8),
    xticklabels=True,
    yticklabels=False,
    dendrogram_ratio=(0.1, 0.2),
    cbar_pos=(1, 0.3, 0.03, 0.5)
)


# Rotate x-axis labels
plt.setp(clustermap.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

# Save figure
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/pathway_clustermap_100.png", dpi=300, bbox_inches="tight")  
plt.show()




# Save plasmid names grouped by cluster (FIXED)
cluster_groups = binary_matrix_filtered.reset_index().groupby("Cluster")["Plasmid"].apply(lambda x: ", ".join(x))

# Save to a CSV file
cluster_groups.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_100_cl_26.csv", header=True)

print("Plasmid names per cluster saved to 'plasmid_cluster_mapping_100_cl_26.csv'")

# Load the cluster mapping file
cluster_mapping = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_100_cl_26.csv", 
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

# Filter clusters with more than 10 plasmids
filtered_clusters = cluster_mapping[cluster_mapping["Plasmid_Count"] > 20]

# Save filtered clusters
filtered_clusters.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_filtered_100_cl_26.csv", 
    index=True,  
    header=True,  
    columns=["Plasmid_Count", "Percentage"]
)

print("\nFiltered clusters (more than 20 plasmids) saved to 'plasmid_cluster_filtered_100_cl_26.csv'.")

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
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_cluster_barplot_100_cl_26.png", dpi=300, bbox_inches="tight")
plt.show()

# Reset index so that "Plasmid" is a column instead of an index
binary_matrix_filtered = binary_matrix_filtered.reset_index()

# Extract metagenome name from sequence names
binary_matrix_filtered["Metagenome"] = binary_matrix_filtered["Plasmid"].str.split("_").str[0]

# Load metadata
metadata = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt", 
    sep="\t", 
    names=["Run", "biome_genus"],  
    header=0  # Adjust if needed
)

binary_matrix_filtered = binary_matrix_filtered.merge(
    metadata, left_on="Metagenome", right_on="Run", how="left", suffixes=("", "_drop")
)

# Drop unnecessary columns
binary_matrix_filtered = binary_matrix_filtered.drop(columns=["Run", "biome_genus_drop"], errors="ignore")


# Verify the change
print(binary_matrix_filtered.head())

# Create color mapping for hosts
unique_hosts = binary_matrix_filtered["biome_genus"].dropna().unique()
num_hosts = len(unique_hosts)
palette = sns.color_palette("hls", n_colors=num_hosts)
host_palette = dict(zip(unique_hosts, palette))


binary_matrix_filtered["Host_Color"] = binary_matrix_filtered["biome_genus"].map(host_palette)

# Select only numeric columns for clustering
# Remove columns that shouldn't be in the heatmap
columns_to_exclude = ["Cluster", "Host_Color", "biome_genus", "Metagenome", "Plasmid"]
numeric_matrix = binary_matrix_filtered.drop(columns=columns_to_exclude, errors="ignore")
numeric_matrix = numeric_matrix.select_dtypes(include=["number"]).copy()


norm = Normalize(vmin=0, vmax=np.max(numeric_matrix))

# Create clustermap with host colors on the left
row_colors = binary_matrix_filtered["Host_Color"]
row_colors.name= "Host"
clustermap = sns.clustermap(
    numeric_matrix,
    cmap="Greys",
    norm=norm,
    row_cluster=True,
    col_cluster=True,
    figsize=(15, 10),
    xticklabels=True,
    yticklabels=False,
    dendrogram_ratio=(0.2, 0.1),
    cbar_pos=None,
    method="ward",
    #row_colors=row_colors  # LEFT side host bar
)

# Add host legend (bottom-left)
legend_patches = [
    mpatches.Patch(color=color, label=host) for host, color in host_palette.items()
]
clustermap.ax_heatmap.legend(
    handles=legend_patches,
    title="Host Genus",
    bbox_to_anchor=(-0.3, -0.2),
    loc="upper left",
    borderaxespad=0,
    fontsize=8,
    title_fontsize=9,
    ncol=2
)

# Extract clusters
linkage_matrix = clustermap.dendrogram_row.linkage
cluster_cutoff = 26
clusters = fcluster(linkage_matrix, cluster_cutoff, criterion="distance")

# Reorder clusters to match heatmap
row_order = clustermap.dendrogram_row.reordered_ind
clusters_ordered = clusters[row_order]

# Create color bar for clusters (RIGHT side)
unique_clusters = sorted(np.unique(clusters_ordered))
cluster_palette = sns.color_palette("colorblind", n_colors=len(unique_clusters))
cluster_color_map = {c: cluster_palette[i] for i, c in enumerate(unique_clusters)}
cluster_colors_ordered = [cluster_color_map[c] for c in clusters_ordered]

# Draw manual cluster bar
pos_heatmap = clustermap.ax_heatmap.get_position()
cluster_bar_ax = clustermap.fig.add_axes([
    pos_heatmap.x1 + 0.02,
    pos_heatmap.y0,
    0.02,
    pos_heatmap.height
])
cluster_bar_ax.imshow(np.array(cluster_colors_ordered).reshape(-1, 1, 3), aspect='auto', origin='upper')
cluster_bar_ax.set_xticks([])
cluster_bar_ax.set_yticks([])
cluster_bar_ax.set_ylabel("Cluster", fontsize=10)
cluster_bar_ax.yaxis.set_label_position("right")
cluster_bar_ax.yaxis.labelpad = 20


# Add cluster number labels next to bar
for cluster in unique_clusters:
    indices = np.where(clusters_ordered == cluster)[0]
    median_idx = indices[len(indices) // 2]
    y_pos = pos_heatmap.y0 + pos_heatmap.height * (1 - median_idx / len(clusters_ordered))
    clustermap.fig.text(
        cluster_bar_ax.get_position().x1 + 0.005,
        y_pos,
        str(cluster),
        va="center",
        ha="left",
        fontsize=9,
        color=cluster_color_map[cluster],
        fontweight="bold"
    )

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/pathway_clustermap_with_clusters.png", dpi=300, bbox_inches="tight", pad_inches=0.1)
plt.show()

# Save cluster assignments
binary_matrix_filtered["Cluster"] = clusters
binary_matrix_filtered[["Plasmid", "Cluster"]].to_csv(
    "C:/Users/hayat/Downloads/R_files/data/Plasmid_Clusters_pathway_clustermap.csv", index=False
)
print("Cluster assignments saved.")




# Save cluster assignments to CSV ===
output_path = "C:/Users/hayat/Downloads/R_files/data/Plasmid_Clusters_pathway_clustermap.csv"
binary_matrix_filtered[["Plasmid", "Cluster"]].to_csv(output_path, index=False)
print(f"Cluster assignments saved to {output_path}")
