# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 14:11:48 2025

@author: hayat
"""

import pandas as pd
import umap
import hdbscan
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
#import numpy as np
#import umap.plot
# Load data
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/KEGG_pathway_per_plasmid.txt", sep="\t", header=None, names=["Plasmid", "Pathway"])

# Convert Pathway to list
df["Pathway"] = df["Pathway"].apply(lambda x: x.split(","))

# Create binary presence/absence matrix
binary_matrix = df.explode("Pathway").assign(Value=1).pivot(index="Plasmid", columns="Pathway", values="Value").fillna(0)

# Filter Pathway that appear in at least 10 plasmids
filtered_columns = binary_matrix.columns[binary_matrix.sum(axis=0) >= 100]
binary_matrix_filtered = binary_matrix[filtered_columns]

# Standardize the data
scaler = StandardScaler()
binary_matrix_std = scaler.fit_transform(binary_matrix_filtered)

# HDBSCAN clustering before UMAP
hdbscan_model = hdbscan.HDBSCAN(min_cluster_size=8, min_samples=5, cluster_selection_epsilon=1.0)
hdbscan_labels = hdbscan_model.fit_predict(binary_matrix_std)

# Add HDBSCAN cluster labels to the data
binary_matrix_filtered["HDBSCAN_Cluster"] = hdbscan_labels

# UMAP dimensionality reduction
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42, densmap=True)
umap_embedding = umap_model.fit_transform(binary_matrix_std)

# Convert UMAP results into a DataFrame
umap_df = pd.DataFrame(umap_embedding, columns=["UMAP1", "UMAP2"], index=binary_matrix_filtered.index).reset_index()
umap_df["HDBSCAN_Cluster"] = hdbscan_labels

# # Plot UMAP with HDBSCAN clusters
# plt.figure(figsize=(7, 5))
# sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue=umap_df["HDBSCAN_Cluster"].astype(str), palette="Dark2", alpha=0.7)
# plt.title("HDBSCAN Clustering of UMAP Results")
# plt.legend(title="Cluster")
# plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_HDBSCAN_plot_Pathway_std.png", dpi=300, bbox_inches='tight')
# plt.show()

# DBSCAN clustering on UMAP
dbscan_model = DBSCAN(eps=1, min_samples=10)
umap_df["DBSCAN_Cluster"] = dbscan_model.fit_predict(umap_df[["UMAP1", "UMAP2"]])

# Compute cluster centers
cluster_centers = umap_df.groupby("DBSCAN_Cluster")[["UMAP1", "UMAP2"]].mean()

# Get unique clusters
clusters = sorted(umap_df["DBSCAN_Cluster"].unique())

# Assign black for outliers (-1), then generate distinct colors for remaining groups
num_clusters = len(clusters) - ("-1" in clusters)  # Exclude outliers when picking colors
color_list = sns.color_palette("turbo", n_colors=num_clusters)  # Try "hsv", "tab20" or "Spectral" if needed

# Create a custom palette dictionary
custom_palette = {"-1": "black"}  # Outlier as black
for i, cluster in enumerate([c for c in clusters if c != -1]):  # Assign colors to other clusters
    custom_palette[str(cluster)] = color_list[i]

# Convert cluster labels to strings (important for seaborn)
umap_df["DBSCAN_Cluster"] = umap_df["DBSCAN_Cluster"].astype(str)

# Plot UMAP with custom colors
plt.figure(figsize=(7, 5))
sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="DBSCAN_Cluster", 
                palette=custom_palette, alpha=0.7, s=10)

# Add cluster labels at centroids
for cluster, (x, y) in cluster_centers.iterrows():
    plt.text(x, y, str(cluster), fontsize=8, ha='center', va='center', fontweight='bold', color='black')

plt.title("DBSCAN Clustering of UMAP Results according to Pathway")

# Hide legend if too large
plt.legend([], [], frameon=False)

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_DBSCAN_plot_pathway_std_py_100.png", dpi=300, bbox_inches='tight')
plt.show()




# Get list of plasmids per cluster
plasmid_list = umap_df.groupby("DBSCAN_Cluster")["Plasmid"].apply(list).to_dict()

# Save plasmid list to text file
with open("C:/Users/hayat/Downloads/R_files/data/plasmids_per_cluster_Pathway_std_100.txt", "w") as f:
    for cluster, plasmids in plasmid_list.items():
        f.write(f"Cluster {cluster}:\n")
        f.write("\n".join(plasmids) + "\n\n")

# Load host information
host_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt", sep="\t")

# Extract metagenome ID from plasmid names
df["Metagenome_ID"] = df["Plasmid"].apply(lambda x: x.split("_")[0])

# Merge with host information
df = df.merge(host_df, left_on="Metagenome_ID", right_on="Run", how="left")
binary_matrix_filtered["Host"] = df.set_index("Plasmid")["biome_genus"]

# Generate a color palette with distinct colors
num_hosts = df["biome_genus"].nunique()
palette = sns.color_palette("tab20", n_colors=num_hosts)  # Use 'husl' if you have more than 20 categories

plt.figure(figsize=(7, 5))
sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue=df["biome_genus"], palette=palette, alpha=0.7)

# Add cluster labels at centroids
for cluster, (x, y) in cluster_centers.iterrows():
    plt.text(x, y, str(cluster), fontsize=8, ha='center', va='center', fontweight='bold', color='black')

plt.title("UMAP Plot Colored by Host Genus")

# Adjust legend to have multiple columns
plt.legend(title="Host Genus", bbox_to_anchor=(1.05, 1), loc="upper left", ncol=3)  # Change ncol to 2 or more as needed

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_Host_Color_pathway_100.png", dpi=300, bbox_inches='tight')
plt.show()


umap_df = umap_df.merge(df[["Plasmid", "biome_genus"]], left_on="Plasmid", right_on="Plasmid", how="left")

umap_df["biome_genus"] = umap_df["biome_genus"].fillna("Aplysina")


g = sns.FacetGrid(umap_df, col="biome_genus", col_wrap=5, height=4, sharex=True, sharey=True)
g.map(sns.scatterplot, "UMAP1", "UMAP2", alpha=0.7)
g.set_titles("{col_name}")
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/UMAP_Facet_Grid_Host_pathway_100.png", dpi=300, bbox_inches='tight')
plt.show()



