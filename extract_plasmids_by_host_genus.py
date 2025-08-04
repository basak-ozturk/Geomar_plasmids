# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 11:05:32 2025

@author: hayat
"""

import pandas as pd

# === Load data ===
mapping_file = "C:/Users/hayat/Downloads/R_files/data/sponge_plasmids_rpkm_submatrix_6945.tsv"  
metadata_file = "C:/Users/hayat/Downloads/R_files/data/GEOLOCATION_HostGenus_Samples.txt" 

# Read mapping matrix
mapping_df = pd.read_csv(mapping_file, sep='\t', index_col=0)

# Read metadata (Run = metagenome name)
meta_df = pd.read_csv(metadata_file, sep='\t')

# Target biome genera
target_genera = ["Aplysina"]

# Process each genus
for genus in target_genera:
    # Step 1: Get metagenomes from this genus
    genus_runs = meta_df[meta_df["biome_genus"] == genus]["Run"].tolist()
    genus_runs = [run for run in genus_runs if run in mapping_df.columns]  # ensure they're in mapping

    if not genus_runs:
        print(f"No metagenomes found for genus {genus}. Skipping.")
        continue

    # Step 2: Extract submatrix
    sub_df = mapping_df[genus_runs].copy()

    # Step 3: Filter plasmids with any RPKM >= 1
    mask = (sub_df >= 1).any(axis=1)
    filtered_df = sub_df[mask]

    # Step 4: Build TSV output
    output_rows = []

    for plasmid, row in filtered_df.iterrows():
        present_in = row[row >= 1]
        metagenomes = ";".join(present_in.index)
        total_rpkm = present_in.sum()
        output_rows.append([plasmid, metagenomes, total_rpkm])

    output_df = pd.DataFrame(output_rows, columns=["Plasmid", "Metagenomes", "Total_RPKM"])

    # Step 5: Save
    output_filename = f"C:/Users/hayat/Downloads/R_files/data/plasmids_in_{genus}.tsv"
    #output_df.to_csv(output_filename, sep="\t", index=False)
    print(f"Saved: {output_filename}")
