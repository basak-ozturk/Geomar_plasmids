# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 16:46:01 2025

@author: hayat
"""

import pandas as pd

# Load your TSV file
input_file = "C:/Users/hayat/Downloads/R_files/data/widespread_plasmid_ani_large_clusters.tsv"
df = pd.read_csv(input_file, sep="\t")

# Define cluster colors
cluster_colors = {
    2: "#e41a1c",     # red
    11: "#377eb8",    # blue
    12: "#4daf4a",    # green
    13: "#984ea3",    # purple
    21: "#ff7f00"     # orange
}

# Create header for iTOL annotation
header = [
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tANI Clusters",
    "COLOR\t#ff0000",  # Default strip color
    "",
    "LEGEND_TITLE\tANI Cluster",
    "LEGEND_SHAPES\t" + "\t".join(["1"] * len(cluster_colors)),
    "LEGEND_COLORS\t" + "\t".join(cluster_colors[c] for c in sorted(cluster_colors)),
    "LEGEND_LABELS\t" + "\t".join(f"Cluster {c}" for c in sorted(cluster_colors)),
    "",
    "DATA"
]

# Generate annotation lines
data_lines = [
    f"{row['Plasmid']}\t{cluster_colors[row['Cluster']]}\tCluster {row['Cluster']}"
    for _, row in df.iterrows()
]

# Combine header and data
itol_lines = header + data_lines

# Save to file
output_file = "C:/Users/hayat/Downloads/R_files/data/itol_ani_clusters.txt"
with open(output_file, "w") as f:
    f.write("\n".join(itol_lines))

print(f"iTOL annotation file saved to: {output_file}")
