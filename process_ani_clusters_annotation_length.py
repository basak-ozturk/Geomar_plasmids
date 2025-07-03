# -*- coding: utf-8 -*-
"""
Created on 2025-07-01

Full pipeline for:
- Extracting large plasmid clusters (≥5 members)
- Filtering eggNOG annotations per cluster
- Extracting FASTA sequences per cluster
- Calculating length stats per cluster
"""

import pandas as pd
import os
from Bio import SeqIO
import numpy as np

# === File paths ===
original_cluster_file = "C:/Users/hayat/Downloads/R_files/data/widespread_plasmid_clusters_cutoff_75.tsv"  # full cluster assignments
large_cluster_file = "C:/Users/hayat/Downloads/R_files/data/widespread_plasmid_ani_large_clusters.tsv"
annotations_file = "C:/Users/hayat/Downloads/R_files/data/filtered_foreground_eggnog_annotations.tsv"
fasta_file = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid.fasta"
cluster_annotation_folder = "C:/Users/hayat/Downloads/R_files/data/ani_cluster_annotations"

# Create output folder if missing
os.makedirs(cluster_annotation_folder, exist_ok=True)

# === Step 1: Load original clusters and extract large clusters (≥5 plasmids) ===
print("Loading original clusters...")
clusters = pd.read_csv(original_cluster_file, sep="\t")

print("Counting plasmids per cluster...")
cluster_counts = clusters['Cluster'].value_counts()

print("Filtering clusters with >= 5 plasmids...")
large_clusters = cluster_counts[cluster_counts >= 5].index

large_clusters_df = clusters[clusters['Cluster'].isin(large_clusters)]

print(f"Saving filtered large clusters to {large_cluster_file} ...")
large_clusters_df.to_csv(large_cluster_file, sep="\t", index=False)

# === Step 2: Filter annotations and save per-cluster annotation files ===
print("Loading eggNOG annotations...")
annotations = pd.read_csv(annotations_file, sep="\t")

annotations['Plasmid'] = annotations['#query'].apply(lambda x: "_".join(x.split("_")[:-1]))

print("Filtering annotations for plasmids in large clusters...")
annotated_subset = annotations[annotations['Plasmid'].isin(large_clusters_df['Plasmid'])]

print("Merging annotations with cluster info...")
annotated_with_cluster = annotated_subset.merge(large_clusters_df, on='Plasmid')

print("Saving one annotation file per cluster...")
for cluster_id, group in annotated_with_cluster.groupby("Cluster"):
    output_path = os.path.join(cluster_annotation_folder, f"cluster_{cluster_id}.tsv")
    group.to_csv(output_path, sep="\t", index=False)

# === Step 3: Load FASTA sequences and index by plasmid name ===
print("Indexing FASTA records...")
fasta_index = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# === Step 4: Extract FASTA sequences and calculate length stats per cluster ===
stats_list = []

print("Extracting FASTA sequences and calculating stats per cluster...")
for cluster_id, group in large_clusters_df.groupby("Cluster"):
    plasmids = group['Plasmid'].unique()
    
    # Get all records for this cluster
    records = [fasta_index[p] for p in plasmids if p in fasta_index]
    
    # Save FASTA file per cluster
    output_fasta = os.path.join(cluster_annotation_folder, f"cluster_{cluster_id}.fasta")
    if records:
        SeqIO.write(records, output_fasta, "fasta")
    else:
        print(f"[Warning] No sequences found in FASTA for cluster {cluster_id}")
    
    # Calculate length stats
    lengths = [len(rec.seq) for rec in records]
    if lengths:
        stats_list.append({
            "Cluster": cluster_id,
            "Num_Plasmids": len(lengths),
            "Total_Length": sum(lengths),
            "Mean_Length": round(np.mean(lengths), 2),
            "Median_Length": round(np.median(lengths), 2),
            "Min_Length": min(lengths),
            "Max_Length": max(lengths)
        })
    else:
        stats_list.append({
            "Cluster": cluster_id,
            "Num_Proteins": 0,
            "Total_Length": 0,
            "Mean_Length": None,
            "Median_Length": None,
            "Min_Length": None,
            "Max_Length": None
        })

# === Step 5: Save length stats summary ===
stats_df = pd.DataFrame(stats_list)
stats_df.sort_values("Cluster", inplace=True)

length_stats_file = os.path.join(cluster_annotation_folder, "cluster_length_stats.tsv")
print(f"Saving length stats to {length_stats_file} ...")
stats_df.to_csv(length_stats_file, sep="\t", index=False)

print("Done!")
