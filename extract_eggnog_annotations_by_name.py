# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 11:24:44 2025

@author: hayat
"""

import pandas as pd

from Bio import SeqIO

# Load plasmid names (no suffix)
with open("C:/Users/hayat/Downloads/R_files/data/itol_node_I187.txt") as f:
    plasmids = set(line.strip() for line in f)

# Read eggNOG annotation TSV
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/eggnog_top_widespread.tsv", sep="\t")

# Extract plasmid part before the last underscore
df["Plasmid"] = df["#query"].str.extract(r"(^.+_\d+ctg)")

# Filter for those in the plasmid list
filtered_df = df[df["Plasmid"].isin(plasmids)]

# Optional: drop the helper column
filtered_df = filtered_df.drop(columns=["Plasmid"])

# Save or inspect
#filtered_df.to_csv("C:/Users/hayat/Downloads/R_files/data/filtered_annotations_node_I187.tsv", sep="\t", index=False)


# Input and output FASTA files
input_fasta = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid.fasta"
output_fasta = "C:/Users/hayat/Downloads/R_files/data/itol_node_I187.fasta"

# Filter and write
with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Match by exact name before any underscores or as-is
        record_id = record.id.split()[0]  # only the first word in header
        if record_id in plasmids:
            SeqIO.write(record, out_f, "fasta")