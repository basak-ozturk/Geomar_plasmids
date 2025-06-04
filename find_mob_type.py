# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 09:00:28 2025

@author: hayat
"""

#import os
import re
import csv
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Directory containing your MOB result files
input_dir = Path(r"C:\Users\hayat\Downloads\R_files\data\conjscan\hmmer_results")

# Pattern to match MOB result files
mob_pattern = re.compile(r"T4SS_MOB(\w+)\.res_hmm_extract")

# Prepare output list
mob_hits = []

# Scan each MOB file
for file in input_dir.glob("T4SS_MOB*.res_hmm_extract"):
    match = mob_pattern.search(file.name)
    if not match:
        continue
    mob_type = match.group(1)

    with open(file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if parts:
                gene_id = parts[0]
                gene_base = "_".join(gene_id.split("_")[:-1])
                mob_hits.append((gene_base, f"MOB{mob_type}"))

# Remove duplicates and sort
mob_hits = sorted(set(mob_hits))

# Save to CSV
output_file = input_dir / "mob_relaxase_hits.csv"
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Plasmid", "MOB type"])
    writer.writerows(mob_hits)

print(f"Saved results to: {output_file}")


# Load the CSV
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/conjscan/hmmer_results/mob_relaxase_hits.csv")

# Count unique plasmid names
unique_count = df["Plasmid"].nunique()

print(f"Unique plasmid names: {unique_count}")

# Count each MOB type
mob_counts = df["MOB type"].value_counts().sort_values(ascending=False)

# Plot
plt.figure(figsize=(10, 6))
mob_counts.plot(kind='bar', color='skyblue', edgecolor='black')
plt.title("Frequency of MOB Types", fontsize=14)
plt.xlabel("MOB Type", fontsize=12)
plt.ylabel("Count", fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Extract metagenome ID from the plasmid name
df['metagenome'] = df['Plasmid'].str.extract(r'^([^_]+)_')

# Create a presence/absence matrix
co_occurrence = df.drop_duplicates()[['metagenome', 'MOB type']]
co_matrix = pd.crosstab(co_occurrence['metagenome'], co_occurrence['MOB type'])

plt.figure(figsize=(12, 6))
sns.heatmap(co_matrix, cmap="Blues", cbar=True, linewidths=0.5)
plt.title("MOB Type Presence per Metagenome", fontsize=14)
plt.xlabel("MOB Type")
plt.ylabel("Metagenome")
plt.tight_layout()
plt.show()


sponge_map = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt", sep="\t")
sponge_map.columns = ['Run', 'biome_genus']


# Extract metagenome ID from plasmid name
df['metagenome'] = df['Plasmid'].str.extract(r'^([^_]+)_')

# Merge with sponge genus mapping
df_merged = df.merge(sponge_map, left_on='metagenome', right_on='Run', how='left')

# Drop duplicate MOB hits per metagenome per MOB type (optional, avoids double counting)
unique_hits = df_merged.drop_duplicates(subset=['metagenome', 'MOB type'])

# Group by sponge genus and MOB type
mob_by_genus = unique_hits.groupby(['biome_genus', 'MOB type']).size().unstack(fill_value=0)

# Stacked bar plot
mob_by_genus.plot(kind='bar', stacked=True, figsize=(12, 6), colormap='tab20')

ax = mob_by_genus.plot(
    kind='bar', 
    stacked=True, 
    figsize=(12, 6), 
    colormap='tab20', 
    align='center'  # explicitly center bars on ticks
)

ax.set_title("Distribution of MOB Types per Sponge Genus", fontsize=14)
ax.set_xlabel("Sponge Genus", fontsize=12)
ax.set_ylabel("Number of MOB Hits", fontsize=12)

# Rotate ticks and align them correctly
ax.set_xticks(range(len(mob_by_genus.index)))  # ensure ticks match bar positions
ax.set_xticklabels(mob_by_genus.index, rotation=45, ha='right')

plt.tight_layout()
plt.show()

plt.show()

# Optionally normalize rows (relative MOB composition per genus)
mob_normalized = mob_by_genus.div(mob_by_genus.sum(axis=1), axis=0)

# Create clustered heatmap

g = sns.clustermap(
    mob_normalized,
    cmap="Blues",
    figsize=(10, 8),
    annot=True, fmt=".2f",
    linewidths=0.5
)

# Set axis labels
g.ax_heatmap.set_ylabel("Sponge genus", fontsize=12)
g.ax_heatmap.set_xlabel("MOB type", fontsize=12)

# Add title (set a bit higher to avoid cropping)
plt.suptitle("Proportional Distribution of MOB Types by Sponge Genus", y=1.05, fontsize=14)

# Save with padding to include the title
plt.savefig(
    "C:/Users/hayat/Downloads/R_files/graphs/plasmid_mob_types_all.png",
    dpi=300,
    bbox_inches='tight'
)

plt.show()

### Repeat calculation for top widespread and abundant plasmids

# Load the top plasmid list
top_plasmids = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid_names.csv", header=None)
top_plasmids.columns = ["Plasmid"]

# Filter the MOB dataframe to include only top plasmids
df_mob_top = df[df["Plasmid"].isin(top_plasmids["Plasmid"])]

# Extract metagenome ID from plasmid name (if not already a column)
df_mob_top["metagenome"] = df_mob_top["Plasmid"].str.extract(r"^([^_]+)")

# Load sponge genus mapping
sponge_map = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt", sep="\t")

# Merge to get sponge genus
df_mob_top = df_mob_top.merge(sponge_map, left_on="metagenome", right_on="Run")

# Count MOBs per sponge genus
mob_by_genus_top = df_mob_top.groupby(["biome_genus", "MOB type"]).size().unstack(fill_value=0)

# Normalize by row to get proportions
mob_normalized_top = mob_by_genus_top.div(mob_by_genus_top.sum(axis=1), axis=0)

g = sns.clustermap(
    mob_normalized_top,
    cmap="Blues",
    figsize=(10, 8),
    annot=True, fmt=".2f",
    linewidths=0.5
)

g.ax_heatmap.set_ylabel("Sponge genus", fontsize=12)
g.ax_heatmap.set_xlabel("MOB type", fontsize=12)

plt.suptitle("MOB Types in Top Abundant & Widespread Plasmids by Sponge Genus", y=1.05, fontsize=14)

plt.savefig(
    "C:/Users/hayat/Downloads/R_files/graphs/mob_types_top_plasmids.png",
    dpi=300,
    bbox_inches='tight'
)
plt.show()
