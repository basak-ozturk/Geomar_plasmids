import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.colors as mcolors

df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/COGs_per_plasmid_cleaned.txt", 
                 sep=" ",  # Changed from \t to space
                 header=None, 
                 names=["Plasmid", "COG"])

# Clean COG lists after splitting
df["COG"] = df["COG"].apply(
    lambda x: [cog.strip() for cog in str(x).split(",") if cog.strip()] 
    if pd.notna(x) else []
)


# Load COG categories
cog_categories = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/COG_categories.txt", sep="\t", header=None, names=["Category", "Description"])

# Filter out COGs belonging to categories X and S
excluded_categories = {"X", "S"}
cog_categories_filtered = cog_categories[~cog_categories["Category"].isin(excluded_categories)]
filtered_categories = set(cog_categories_filtered["Category"])


def filter_cogs(cog_list):
    valid_cogs = [cog for cog in cog_list 
                 if (cog in filtered_categories) and (len(cog) == 1)]
    #print(f"Filtered: {list(set(cog_list) - set(valid_cogs))}")  # Debug
    return valid_cogs

df["COG"] = df["COG"].apply(filter_cogs)

# Remove rows with empty COG lists
df = df[df["COG"].apply(len) > 0]


# Convert COG column back to comma-separated strings for easier processing
df["COG"] = df["COG"].apply(lambda x: ",".join(x))

# Create binary presence/absence matrix
df["COG"] = df["COG"].apply(lambda x: x.split(","))
binary_matrix = df.explode("COG").assign(Value=1).pivot(index="Plasmid", columns="COG", values="Value").fillna(0)

# Filter COGs that appear in at least 100 plasmids
filtered_columns = binary_matrix.columns[binary_matrix.sum(axis=0) >= 100]
binary_matrix_filtered = binary_matrix[filtered_columns]

# Map COG categories to their descriptions for better labeling
cog_dict = dict(zip(cog_categories["Category"], cog_categories["Description"]))
binary_matrix_filtered = binary_matrix_filtered.rename(columns=lambda x: f"{cog_dict.get(x, 'Unknown')} ({x})")

# Ensure "Cluster" is removed if it exists (to avoid conflicts)
if "Cluster" in binary_matrix_filtered.columns:
    binary_matrix_filtered = binary_matrix_filtered.drop(columns=["Cluster"])

# Convert to boolean for memory efficiency (optional but recommended)
binary_matrix_filtered = binary_matrix_filtered.astype(bool)

# Compute distance matrix 
distance_matrix = pdist(binary_matrix_filtered, metric="euclidean")

# Perform hierarchical clustering using average linkage (better for Jaccard)
linkage_matrix = linkage(distance_matrix, method="ward")

# Dynamically calculate a cutoff threshold based on the linkage matrix
cluster_cutoff = 0.7 * linkage_matrix[-10, 2]  # Adjust multiplier (e.g., 0.7) based on dataset

# Assign cluster labels based on the cutoff threshold
cluster_labels = fcluster(linkage_matrix, t=cluster_cutoff, criterion="distance")
cluster_labels_series = pd.Series(cluster_labels, index=binary_matrix_filtered.index, name="Cluster")

# Add cluster labels to the dataframe for downstream analysis (if needed)
binary_matrix_with_clusters = binary_matrix_filtered.copy()
binary_matrix_with_clusters["Cluster"] = cluster_labels_series

# Plot dendrogram with the computed linkage matrix
plt.figure(figsize=(12, 6))
dendrogram(linkage_matrix, labels=binary_matrix_filtered.index.tolist(), color_threshold=cluster_cutoff)
plt.axhline(y=cluster_cutoff, color='r', linestyle='--')
plt.title("Hierarchical Clustering Dendrogram Based on COGs")
#psavlt.xlabel(f"Plasmids (n={len(binary_matrix_filtered)})")
plt.ylabel("Distance")
plt.xticks([])
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/COG_dendrogram_100.png", dpi=300, bbox_inches="tight")
plt.show()

# Clustermap visualization using the precomputed linkage matrix
clustermap = sns.clustermap(
    binary_matrix_filtered,
    row_linkage=linkage_matrix,
    col_cluster=True,
    cmap="Greys",
    figsize=(12, 8),
    xticklabels=1,
    yticklabels=False,
    dendrogram_ratio=(0.1, 0.2),  # Adjust dendrogram size as needed
    cbar=False,
    cbar_pos=None       # Prevent space being reserved for it
)

#plt.show()

# Rotate x-axis labels for better readability in clustermap
#plt.setp(clustermap.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

# Save clustermap figure
clustermap.savefig("C:/Users/hayat/Downloads/R_files/graphs/COG_clustermap_100.png", dpi=300, bbox_inches="tight")
plt.show()

# Reset index to get "Plasmid" as a column
binary_matrix_filtered = binary_matrix_filtered.reset_index()

# Rebuild cluster labels Series with the new index
cluster_labels_series = pd.Series(cluster_labels, index=binary_matrix_filtered["Plasmid"], name="Cluster")

# Add cluster labels as a new column (matches on Plasmid)
binary_matrix_filtered["Cluster"] = binary_matrix_filtered["Plasmid"].map(cluster_labels_series)

# Group by cluster and join plasmid names
cluster_groups = binary_matrix_filtered.groupby("Cluster")["Plasmid"].apply(lambda x: ", ".join(x))

# Save to CSV
cluster_groups.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_COG_100.csv", header=True)

print(cluster_groups.head())
print(binary_matrix_filtered[["Plasmid", "Cluster"]].head())

print("Plasmid names per cluster saved to 'plasmid_cluster_mapping_COG_100.csv'")

# Load the cluster mapping file
cluster_mapping = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_COG_100.csv", 
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
filtered_clusters = cluster_mapping[cluster_mapping["Plasmid_Count"] > 1]

# Save filtered clusters
filtered_clusters.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_filtered_COG_100.csv", 
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
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_cluster_barplot_COG_100.png", dpi=300, bbox_inches="tight")
plt.show()

# Reset index so Plasmid is a column
binary_matrix_filtered = binary_matrix_filtered.reset_index()

# Extract Metagenome from Plasmid
binary_matrix_filtered["Metagenome"] = binary_matrix_filtered["Plasmid"].str.split("_").str[0]

# Load metadata
metadata = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt",
    sep="\t",
    names=["Run", "biome_genus"],
    header=0
)


# Merge and clean up
binary_matrix_filtered = binary_matrix_filtered.merge(
    metadata, left_on="Metagenome", right_on="Run", how="left", suffixes=("", "_drop")
)
binary_matrix_filtered = binary_matrix_filtered.drop(columns=["Run", "biome_genus_drop"], errors="ignore")
binary_matrix_filtered["biome_genus"] = binary_matrix_filtered["biome_genus"].fillna("Aplysina")

# Define host colors
host_palette = {
    host: color for host, color in zip(
        binary_matrix_filtered["biome_genus"].unique(),
        sns.color_palette("tab10", n_colors=binary_matrix_filtered["biome_genus"].nunique())
    )
}
row_colors = binary_matrix_filtered["biome_genus"].map(host_palette)
row_colors.name = ""  # Suppress label in row_colors legend

# Define numeric matrix
numeric_matrix = binary_matrix_filtered.drop(
    columns=["index", "Plasmid", "Metagenome", "biome_genus", "Cluster"],
    errors="ignore"
)
# Extract linkage matrix and form flat clusters
# === Prepare matrix ===

# Ensure clusters are defined before clustermap
linkage_matrix = sns.clustermap(
    numeric_matrix,
    row_cluster=True,
    col_cluster=True,
    cmap="gray_r",
    figsize=(12, 14),
    xticklabels=True,
    yticklabels=False,  # Disable default labels
    row_colors=row_colors,
    cbar_pos=None,
    method="ward",
    dendrogram_ratio=(0.4, 0.1)
).dendrogram_row.linkage

# Add cluster labels before plotting
max_d = 20
clusters = fcluster(linkage_matrix, max_d, criterion="distance")
binary_matrix_filtered["Cluster"] = clusters

# Use bright distinct cluster colors
unique_clusters = sorted(binary_matrix_filtered["Cluster"].unique())
cluster_palette = {
    cluster: color for cluster, color in zip(
        unique_clusters,
        sns.color_palette("tab20", n_colors=len(unique_clusters))
    )
}



# === Generate clustermap ===
clustermap = sns.clustermap(
    numeric_matrix,
    row_cluster=True,
    col_cluster=True,
    cmap="gray_r",
    figsize=(12, 14),
    xticklabels=True,
    yticklabels=False,
    row_colors=row_colors,
    cbar_pos=None,
    method="ward",
    dendrogram_ratio=(0.4, 0.1)
)

# === Centered, unique cluster labels on y-axis ===

# === Step 1: Prepare host colors (unordered) ===
host_color_list = binary_matrix_filtered["biome_genus"].map(host_palette).tolist()

# === Step 2: Initial clustermap to get row order ===
clustermap = sns.clustermap(
    numeric_matrix,
    row_cluster=True,
    col_cluster=True,
    cmap="gray_r",
    figsize=(12, 14),
    xticklabels=True,
    yticklabels=False,
    row_colors=host_color_list,  # <- host on the left
    cbar_pos=None,
    method="ward",
    dendrogram_ratio=(0.4, 0.1)
)
fig = clustermap.fig
ax = clustermap.ax_heatmap
pos = ax.get_position()
x0, y0, width, height = pos.x0, pos.y0, pos.width, pos.height


# === Step 3: Reorder indices and colors ===
reordered_indices = clustermap.dendrogram_row.reordered_ind
cluster_labels_ordered = binary_matrix_filtered.iloc[reordered_indices]["Cluster"].values
host_color_list_ordered = binary_matrix_filtered.iloc[reordered_indices]["biome_genus"].map(host_palette).tolist()
cluster_color_list_ordered = [cluster_palette[int(c)] for c in cluster_labels_ordered]

# === Step 4: Add cluster colorbar manually to the right (corrected version) ===
rgb_colors = [mcolors.to_rgb(cluster_palette[int(c)]) for c in cluster_labels_ordered]

cluster_bar_ax = fig.add_axes([x0 + width + 0.03, y0, 0.02, height])  # placed more to the right
cluster_bar_ax.imshow(
    np.array(rgb_colors).reshape(-1, 1, 3),
    aspect='auto',
    origin='upper'
)
cluster_bar_ax.set_xticks([])
cluster_bar_ax.set_yticks([])
cluster_bar_ax.set_ylabel("Clusters", fontsize=10)
cluster_bar_ax.yaxis.set_label_position("right")

# === Step 5: Center cluster number labels on y-axis ===
unique_clusters, cluster_counts = np.unique(cluster_labels_ordered, return_counts=True)
label_positions = []
for cluster in unique_clusters:
    indices = np.where(cluster_labels_ordered == cluster)[0]
    label_positions.append(indices[len(indices) // 2])  # median position


# === Step 5: Centered labels — moved further right to avoid overlap ===
clustermap.ax_heatmap.set_yticks(label_positions)
clustermap.ax_heatmap.set_yticklabels(unique_clusters, fontsize=9, va="center", ha="left")
for pos, cluster in zip(label_positions, unique_clusters):
    label = clustermap.ax_heatmap.get_yticklabels()[list(unique_clusters).index(cluster)]
    label.set_color(cluster_palette[int(cluster)])
    label.set_fontweight("bold")
    label.set_x(1.02)  # Shift to the right

# Cleanup
clustermap.ax_heatmap.tick_params(axis='y', length=0)
clustermap.ax_heatmap.spines['left'].set_visible(False)

# === Step 6: Legends ===
host_patches = [mpatches.Patch(color=v, label=k) for k, v in host_palette.items()]
cluster_patches = [mpatches.Patch(color=v, label=f"Cluster {k}") for k, v in cluster_palette.items()]

clustermap.ax_heatmap.legend(
    handles=host_patches,
    title="Host Genus",
    loc="upper left",
    bbox_to_anchor=(1.2, 1),
    borderaxespad=0,
    ncol=2
)

# === Step 7: Save & Show ===
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/COG_clustermap_with_hosts_right_clusterbar.png",
            dpi=300, bbox_inches="tight")
plt.show()

# Extract linkage matrix from the clustermap
linkage_matrix = clustermap.dendrogram_row.linkage

# Define the number of clusters or distance threshold (adjust as needed)
max_d = 20  # Adjust based on the dendrogram structure
clusters = fcluster(linkage_matrix, max_d, criterion="distance")

# Add cluster assignments to the dataframe
binary_matrix_filtered["Cluster"] = clusters

# Save to CSV, grouping by cluster
output_path = "C:/Users/hayat/Downloads/R_files/data/Plasmid_Clusters_COG_clustermap.csv"
binary_matrix_filtered[["Plasmid", "Cluster"]].to_csv(output_path, index=False)

print(f"Cluster assignments saved to {output_path}")