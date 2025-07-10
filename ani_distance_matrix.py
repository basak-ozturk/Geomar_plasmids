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
from scipy.cluster.hierarchy import to_tree
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
#dist_matrix.to_csv("C:/Users/hayat/Downloads/R_files/data/ani_distance_matrix.tsv", sep="\t")

# Plot clustered heatmap without tick labels
sns.set(style="white")
plt.figure(figsize=(10, 8))
g = sns.clustermap(dist_matrix, cmap="magma", xticklabels=False, yticklabels=False)
plt.suptitle("Clustered ANI Distance Matrix of the Top Abundant and Widespread Plasmids (100 - ANI%)", y=1.02)
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/clustered_ani_heatmap_no_labels.png", dpi=300, bbox_inches='tight')
plt.show()

# Convert to condensed format for linkage
condensed_dist = squareform(dist_matrix.values)
linkage_matrix = linkage(condensed_dist, method='average')

# Plot dendrogram
plt.figure(figsize=(10, 6))
dendrogram(linkage_matrix, labels=dist_matrix.index.tolist(), leaf_rotation=90)
plt.title("Plasmid Dendrogram (ANI Distance, UPGMA)")
plt.xticks([], [])  # Removes both ticks and labels
plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_tree_ani_dendrogram.png", dpi=300)
plt.show()


# ---- Recursive function to convert to Newick ----
def get_newick(node, labels, newick="", parentdist=0.0):
    if node.is_leaf():
        return f"{labels[node.id]}:{parentdist - node.dist:.6f}{newick}"
    else:
        left = get_newick(node.left, labels, "", node.dist)
        right = get_newick(node.right, labels, "", node.dist)
        return f"({left},{right}):{parentdist - node.dist:.6f}{newick}"

# ---- Build tree and generate Newick string ----
tree = to_tree(linkage_matrix, rd=False)
newick_str = get_newick(tree, dist_matrix.index.tolist()) + ";"

# ---- Save to file ----
with open("C:/Users/hayat/Downloads/R_files/data/plasmid_ani_tree.nwk", "w") as f:
    f.write(newick_str)
    
# expected = set(plasmids)  
# present = set(dist_matrix.index)  

# missing = expected - present
# print("Missing plasmids:", missing)


# # Check for duplicate rows
# rows = dist_matrix.values
# unique_rows = np.unique(rows, axis=0)
# print(f"Unique rows: {unique_rows.shape[0]}, Total rows: {rows.shape[0]}")

for cutoff in [90, 85, 80, 75, 70, 65, 60, 55, 50, 25, 10]:
    labels = fcluster(linkage_matrix, cutoff, criterion='distance')
    n_clusters = len(set(labels))
    print(f"Cutoff {cutoff}: {n_clusters} clusters")



threshold = 75.0

# Assign cluster labels based on threshold
cluster_labels = fcluster(linkage_matrix, threshold, criterion='distance')

# Map cluster numbers to plasmid names
cluster_df = pd.DataFrame({
    'Plasmid': dist_matrix.index,
    'Cluster': cluster_labels
})

# Save cluster assignments
#cluster_df.to_csv("C:/Users/hayat/Downloads/R_files/data/widespread_plasmid_clusters_cutoff_75.tsv", sep="\t", index=False)

# Print clusters (optional)
for cluster_id in sorted(cluster_df['Cluster'].unique()):
    members = cluster_df[cluster_df['Cluster'] == cluster_id]['Plasmid'].tolist()
    print(f"Cluster {cluster_id} ({len(members)} plasmids):")
    print(", ".join(members))
    print("-" * 40)
    
#Narrow matrix down to the largest clusters   

cluster_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/widespread_plasmid_ani_large_clusters.tsv", sep="\t")

selected_plasmids = cluster_df["Plasmid"].tolist()

sub_matrix = dist_matrix.loc[selected_plasmids, selected_plasmids]

sns.set(style="white")
plt.figure(figsize=(10, 8))
g = sns.clustermap(sub_matrix, cmap="magma", xticklabels=False, yticklabels=False)
plt.suptitle("Clustered ANI Distance Matrix (Selected Plasmids)", y=1.02)
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/clustered_ani_heatmap_sublist.png", dpi=300, bbox_inches='tight')

plt.show()

condensed = squareform(sub_matrix.values)
linkage_matrix = linkage(condensed, method='average')

plt.figure(figsize=(12, 6))
dendrogram(linkage_matrix, labels=sub_matrix.index, leaf_rotation=90)
plt.xticks([], [])  # Removes both ticks and labels
plt.title("Dendrogram of Selected Plasmids")
plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_tree_dendrogram_sublist.png", dpi=300)

plt.show()