# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 08:45:08 2025

@author: hayat
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import to_tree
#from io import StringIO

# ---- Load Mash distances ----
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_distances.tab", sep="\t", header=None)
df.columns = ["seq1", "seq2", "distance", "p", "shared"]

# ---- Function to clean names ----
def clean_name(name):
    return name.split("/")[-1].replace(".fasta", "")

# ---- Clean the plasmid names in the dataframe ----
df["seq1"] = df["seq1"].apply(clean_name)
df["seq2"] = df["seq2"].apply(clean_name)

# ---- Create sorted list of unique cleaned names ----
plasmids = sorted(set(df['seq1']).union(df['seq2']))

# ---- Initialize distance matrix ----
distance_matrix = pd.DataFrame(np.ones((len(plasmids), len(plasmids))),
                               index=plasmids, columns=plasmids)

# ---- Fill in distances ----
for _, row in df.iterrows():
    distance_matrix.loc[row['seq1'], row['seq2']] = row['distance']
    distance_matrix.loc[row['seq2'], row['seq1']] = row['distance']  # symmetric

# ---- Convert to condensed format ----
condensed = squareform(distance_matrix.values)

# ---- Perform clustering (UPGMA) ----
linkage_matrix = linkage(condensed, method='average')  # UPGMA = average linkage

# ---- Recursive function to convert tree to Newick ----
def get_newick(node, labels, newick="", parentdist=0.0):
    if node.is_leaf():
        return f"{labels[node.id]}:{parentdist - node.dist:.6f}{newick}"
    else:
        left = get_newick(node.left, labels, "", node.dist)
        right = get_newick(node.right, labels, "", node.dist)
        return f"({left},{right}):{parentdist - node.dist:.6f}{newick}"

# ---- Convert to Newick ----
tree = to_tree(linkage_matrix, rd=False)
newick_str = get_newick(tree, plasmids) + ";"

# ---- Save to file ----
with open("C:/Users/hayat/Downloads/R_files/data/plasmid_mash_tree.nwk", "w") as f:
    f.write(newick_str)

# ---- Plot dendrogram ----
plt.figure(figsize=(12, 6))
dendrogram(linkage_matrix, labels=distance_matrix.index, leaf_rotation=90, leaf_font_size=6)
plt.title("Plasmid Dendrogram (Mash Distance, UPGMA)")
plt.xticks([], []) 
plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_dendrogram_mash.png", dpi=300)
plt.show()
