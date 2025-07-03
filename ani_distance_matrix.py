# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 14:15:40 2025

@author: hayat
"""

import pandas as pd
#import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import fcluster

# Load FastANI output
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/fastani_output.txt", sep="\t", header=None)
df.columns = ["query", "ref", "ani", "shared_frags", "total_frags"]

# Clean filenames (remove paths and extensions)
df["query"] = df["query"].apply(lambda x: x.split("/")[-1].replace(".fasta", ""))
df["ref"] = df["ref"].apply(lambda x: x.split("/")[-1].replace(".fasta", ""))

# All unique plasmids
plasmids = sorted(set(df["query"]).union(df["ref"]))
dist_matrix = pd.DataFrame(100.0, index=plasmids, columns=plasmids)

# Fill distances (100 - ANI)
for _, row in df.iterrows():
    dist = 100 - row["ani"]
    dist_matrix.loc[row["query"], row["ref"]] = dist
    dist_matrix.loc[row["ref"], row["query"]] = dist

# Save matrix
dist_matrix.to_csv("C:/Users/hayat/Downloads/R_files/data/ani_distance_matrix.tsv", sep="\t")

# Plot clustered heatmap without tick labels
sns.set(style="white")
plt.figure(figsize=(10, 8))
g = sns.clustermap(dist_matrix, cmap="magma", xticklabels=False, yticklabels=False)
plt.suptitle("Clustered ANI Distance Matrix of the Top Abundant and Widespread Plasmids (100 - ANI%)", y=1.02)
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/clustered_ani_heatmap_no_labels.png", dpi=300, bbox_inches='tight')
plt.show()

# Convert to condensed format for linkage
condensed_dist = squareform(dist_matrix.values)
linkage_matrix = linkage(condensed_dist, method='average')

# Plot dendrogram
plt.figure(figsize=(10, 6))
dendrogram(linkage_matrix, labels=dist_matrix.index.tolist(), leaf_rotation=90)
plt.title("Plasmid Tree (Hierarchical Clustering from ANI)")
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_tree_dendrogram.png", dpi=300)
plt.show()

#Define distance threshold (e.g., 20 = plasmids with ≥80% ANI are grouped)
threshold = 75.0

# Assign cluster labels based on threshold
cluster_labels = fcluster(linkage_matrix, threshold, criterion='distance')

# Map cluster numbers to plasmid names
cluster_df = pd.DataFrame({
    'Plasmid': dist_matrix.index,
    'Cluster': cluster_labels
})

# Save cluster assignments
cluster_df.to_csv("C:/Users/hayat/Downloads/R_files/data/widespread_plasmid_clusters_cutoff_75.tsv", sep="\t", index=False)

# Print clusters (optional)
for cluster_id in sorted(cluster_df['Cluster'].unique()):
    members = cluster_df[cluster_df['Cluster'] == cluster_id]['Plasmid'].tolist()
    print(f"Cluster {cluster_id} ({len(members)} plasmids):")
    print(", ".join(members))
    print("-" * 40)