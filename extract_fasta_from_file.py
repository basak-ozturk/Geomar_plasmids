# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 15:11:19 2025

@author: hayat
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

# File paths — change as needed
fasta_file = "C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.fasta"
names_file = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_names.txt"
output_fasta = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid.fasta"
output_csv = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_sequence_lengths.csv"
output_plot = "C:/Users/hayat/Downloads/R_files/graphs/top_abundant_and_widespread_sequence_length_distribution.png"

# Step 1: Read target names
with open(names_file, "r") as f:
    target_names = set(line.strip() for line in f if line.strip())

# Step 2: Parse input FASTA and extract matching sequences
matched_seqs = [seq for seq in SeqIO.parse(fasta_file, "fasta") if seq.id in target_names]

# Step 3: Write matching sequences to FASTA
SeqIO.write(matched_seqs, output_fasta, "fasta")
print(f"Done. Matching sequences written to: {output_fasta}")

# Step 4: Create DataFrame of lengths
df = pd.DataFrame({"Sequence_ID": [seq.id for seq in matched_seqs],
                   "Length": [len(seq) for seq in matched_seqs]})

# Step 5: Write lengths to CSV
df.to_csv(output_csv, index=False)

# Step 6: Plot length distribution
plt.figure(figsize=(8, 6))
plt.hist(df["Length"], bins=30, color="skyblue", edgecolor="black")
plt.title("Length Distribution of Top Abundant And Widespread Plasmids")
plt.xlabel("Sequence Length")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.show()
plt.close()

print(f"Done! Wrote:\n- CSV: {output_csv}\n- Plot: {output_plot}")
