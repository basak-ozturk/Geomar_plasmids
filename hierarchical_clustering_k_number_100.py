import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, Normalize
import numpy as np
import matplotlib.ticker as ticker

# Load K-number descriptions
kegg_definitions = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/eggnog_all_plasmids_k_per_plasmid.tsv", 
    sep="\t", 
    header=None, 
    names=["Cluster", "KNumber", "Score", "Definition"]
)


# Create a dictionary mapping K-numbers to "Definition (K number)"
kegg_dict = {row["KNumber"]: f"{row['Definition']} ({row['KNumber']})" for _, row in kegg_definitions.iterrows()}

# Function to replace K-numbers with their definitions
def replace_k_numbers(k_list):
    return [kegg_dict.get(k, k) for k in k_list]  # Keep K-number if not found in dictionary

# Load KEGG KO per plasmid data
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_ko_per_plasmid.txt", sep="\t", header=None, names=["Plasmid", "KEGG"])

# Convert KEGG column from comma-separated string to list
df["KEGG"] = df["KEGG"].apply(lambda x: x.split(","))

# Replace K-numbers with their definitions
df["KEGG"] = df["KEGG"].apply(replace_k_numbers)

# Convert back to comma-separated string
df["KEGG"] = df["KEGG"].apply(lambda x: ",".join(x))

# Save the updated file (optional)
df.to_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_ko_per_plasmid_with_definitions.txt", sep="\t", index=False, header=False)

print("Replaced K-numbers with descriptions in KEGG KO file.")


# Convert KEGG to list
df["KEGG"] = df["KEGG"].apply(lambda x: x.split(","))

# Create binary presence/absence matrix
binary_matrix = df.explode("KEGG").assign(Value=1).pivot(index="Plasmid", columns="KEGG", values="Value").fillna(0)

# Filter K-numbers that appear in at least 100 plasmids
filtered_columns = binary_matrix.columns[binary_matrix.sum(axis=0) >= 100]
binary_matrix_filtered = binary_matrix[filtered_columns]

# Compute distance matrix (Jaccard distance)
distance_matrix = pdist(binary_matrix_filtered, metric="euclidean")

# Perform hierarchical clustering (change linkage method if needed)
linkage_matrix = linkage(distance_matrix, method="ward")  # Use 'average' for Jaccard

# Define a cut-off threshold for clustering (adjust as needed)
cluster_cutoff = 20  # Adjust based on the dendrogram structure

# Assign cluster labels
cluster_labels = fcluster(linkage_matrix, t=cluster_cutoff, criterion="distance")

# Add cluster labels to dataframe
binary_matrix_filtered["Cluster"] = cluster_labels

# Save cluster assignments
binary_matrix_filtered.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_clusters_hierarchical_KEGG_100.csv")

# Set Seaborn style
sns.set_style("whitegrid")

# Plot dendrogram and save it
plt.figure(figsize=(12, 6))
#dendrogram(linkage_matrix, no_labels=True)  # Show labels if needed
dendrogram(linkage_matrix, no_labels=True, color_threshold=cluster_cutoff)
plt.axhline(y=cluster_cutoff, color='r', linestyle='--')  # Show cut-off threshold
plt.title("Hierarchical Clustering Dendrogram Based on KEGG numbers")
plt.xlabel("Plasmids")
plt.ylabel("Distance")
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/KEGG_dendrogram_100.png", dpi=300, bbox_inches="tight")  # Save as PNG
plt.show()

clustermap = sns.clustermap(
    binary_matrix_filtered.drop(columns=["Cluster"]),
    row_cluster=True,
    col_cluster=True,
    cmap="Greys",
    figsize=(12, 8),
    xticklabels=True,
    yticklabels=False,
    dendrogram_ratio=(0.1, 0.2),  # Reduce dendrogram size
    cbar_pos=(1, 0.3, 0.03, 0.5)  # Move colorbar down
)

# Rotate x-axis labels
plt.setp(clustermap.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

# Save figure
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/KEGG_clustermap_100.png", dpi=300, bbox_inches="tight")  
plt.show()




# Save plasmid names grouped by cluster (FIXED)
cluster_groups = binary_matrix_filtered.reset_index().groupby("Cluster")["Plasmid"].apply(lambda x: ", ".join(x))

# Save to a CSV file
cluster_groups.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_KEGG_100.csv", header=True)

print("Plasmid names per cluster saved to 'plasmid_cluster_mapping_KEGG_100.csv'")

# Load the cluster mapping file
cluster_mapping = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_mapping_KEGG_100.csv", 
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
filtered_clusters = cluster_mapping[cluster_mapping["Plasmid_Count"] > 50]

# Save filtered clusters
filtered_clusters.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/plasmid_cluster_filtered_KEGG_100.csv", 
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
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_cluster_barplot_KEGG_100.png", dpi=300, bbox_inches="tight")
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

# Replace NaN values in 'biome_genus' with "Aplysina"
binary_matrix_filtered["biome_genus"] = binary_matrix_filtered["biome_genus"].fillna("Aplysina")

# Verify the change
print(binary_matrix_filtered.head())


# Create a color mapping for hosts
host_palette = {host: color for host, color in zip(binary_matrix_filtered["biome_genus"].unique(), sns.color_palette("tab10"))}
binary_matrix_filtered["Host_Color"] = binary_matrix_filtered["biome_genus"].map(host_palette)

# Select only numeric columns for clustering
numeric_matrix = binary_matrix_filtered.select_dtypes(include=["number"])

viridis = plt.cm.viridis(np.linspace(0, 1, 256))

# Set the first color (for zeros) to white
viridis[0] = [1, 1, 1, 1]  # RGBA (White)

# Create a new colormap with modified colors
custom_cmap = ListedColormap(viridis)

# Normalize: Ensure zero maps to the first color (white)
norm = Normalize(vmin=0, vmax=np.max(numeric_matrix))

# Define custom tick locations and labels
tick_locations = [i + 0.5 for i in range(14)]  # Centers of the color intervals
tick_labels = [str(i + 1) for i in range(14)]  # Labels from 1 to 14


# Apply the custom colormap in clustermap
clustermap = sns.clustermap(
    numeric_matrix, 
    cmap=custom_cmap, 
    norm=norm,
    row_cluster=True, 
    col_cluster=True, 
    figsize=(12, 8),
    xticklabels=True,
    yticklabels=False,
    row_colors=binary_matrix_filtered["Host_Color"],
    dendrogram_ratio=(0.2, 0.1),
    cbar_pos=(1.05, 0.3, 0.03, 0.5),
    method="ward"
)

# Get the colorbar axis
cbar_ax = clustermap.ax_cbar

# Set the custom tick locations and labels
cbar_ax.set_yticks(tick_locations)
cbar_ax.set_yticklabels(tick_labels)

# Increase contrast for binary data
#clustermap.ax_heatmap.collections[0].set_clim(0, 1)

# Rotate x-axis labels and decrease font size
plt.setp(clustermap.ax_heatmap.get_xticklabels(), ha="right", fontsize=6)  # Adjust font size

# Add legend for host colors
legend_patches = [
    mpatches.Patch(color=color, label=host) for host, color in host_palette.items()
]
clustermap.ax_heatmap.legend(
    handles=legend_patches, 
    title="Host Genus", 
    loc="upper left", 
    bbox_to_anchor=(-0.2, -0.2), 
    borderaxespad=0
)


# Save figure
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/KEGG_clustermap_with_hosts.png", dpi=300, bbox_inches="tight")
plt.show()

# Extract linkage matrix from the clustermap
linkage_matrix = clustermap.dendrogram_row.linkage

# Define the number of clusters or distance threshold (adjust as needed)
max_d = 20  # Adjust based on the dendrogram structure
clusters = fcluster(linkage_matrix, max_d, criterion="distance")

# Add cluster assignments to the dataframe
binary_matrix_filtered["Cluster"] = clusters

# Save to CSV, grouping by cluster
output_path = "C:/Users/hayat/Downloads/R_files/data/Plasmid_Clusters_kegg_clustermap.csv"
binary_matrix_filtered[["Plasmid", "Cluster"]].to_csv(output_path, index=False)

print(f"Cluster assignments saved to {output_path}")



