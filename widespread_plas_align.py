# -*- coding: utf-8 -*-
"""
Created on Wed May 21 14:33:16 2025

@author: hayat
"""

from Bio import SeqIO
import csv
import pandas as pd
import matplotlib.pyplot as plt

# Load list of plasmid names (adjust if header present)
with open("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid_names.csv") as f:
    reader = csv.reader(f)
    next(reader)  # skip header if there is one
    plasmid_names = set(row[0].strip() for row in reader)

# Extract matching sequences
input_fasta = "C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.fasta"
output_fasta = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid.fasta"

count = 0
with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        if record.id in plasmid_names:
            SeqIO.write(record, outfile, "fasta")
            count += 1

print(f"Extracted {count} plasmids.")


# Parse the FASTA file and extract Sequence_ID and Length
lengths = [(record.id, len(record.seq)) for record in SeqIO.parse(output_fasta, "fasta")]

# Create a DataFrame
df = pd.DataFrame(lengths, columns=["Sequence_ID", "Length_bp"])


# Save to CSV
df.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_lengths_widespread_frequent.csv", index=False)

# Print a preview
print(df.head())

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(df["Length_bp"], bins=30, color='skyblue', edgecolor='black')
plt.title("Length Distribution of the Most Frequent and Widespread Plasmid Sequences")
plt.xlabel("Sequence Length (bp)")
plt.ylabel("Frequency")
plt.grid(True)
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/length_distribution_widepsread_frequent_final_plasmids.png", dpi=300, bbox_inches="tight")

plt.show()

