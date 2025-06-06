# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 10:18:23 2025

@author: hayat
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file
cluster_file = r"C:\Users\hayat\Downloads\R_files\data\Plasmid_Clusters_pathway_clustermap.csv"
df = pd.read_csv(cluster_file)

# Group by the Cluster column
cluster_groups = df.groupby("Cluster")["Plasmid"]

# Output directory (optional, can be set to same as input)
output_dir = os.path.dirname(cluster_file)

# Save each group to a separate text file
for cluster_id, plasmids in cluster_groups:
    output_file = os.path.join(output_dir, f"pathway_cluster_{cluster_id}.txt")
    plasmids.to_csv(output_file, index=False, header=False)

print("Files saved successfully.")

conjscan_ids_file = r"C:\Users\hayat\Downloads\R_files\data\unique_CONJscan_hit_ids.txt"

# Load data
clusters_df = pd.read_csv(cluster_file)
with open(conjscan_ids_file, 'r') as f:
    conjscan_ids = [line.strip() for line in f if line.strip()]

# Filter for CONJscan hits
filtered_df = clusters_df[clusters_df["Plasmid"].isin(conjscan_ids)]

# Group by cluster and save each group to a separate file
for cluster_id, group in filtered_df.groupby("Cluster"):
    output_file = os.path.join(output_dir, f"CONJscan_cluster_{cluster_id}.txt")
    group["Plasmid"].to_csv(output_file, index=False, header=False)

print("Clustered CONJscan hit files saved successfully.")


# Get all relevant cluster numbers from the files
cluster_nums = []
for filename in os.listdir(output_dir):
    if filename.startswith("pathway_cluster_") and filename.endswith(".txt"):
        try:
            cluster_num = int(filename.split("_")[-1].split(".")[0])
            cluster_nums.append(cluster_num)
        except ValueError:
            continue

# Sort for consistent order
cluster_nums = sorted(set(cluster_nums))

# Compare overlaps
for cluster in cluster_nums:
    path_pathway = os.path.join(output_dir, f"pathway_cluster_{cluster}.txt")
    path_conjscan = os.path.join(output_dir, f"CONJscan_cluster_{cluster}.txt")

    if os.path.exists(path_conjscan):  # Only compare if CONJscan file exists
        # Load sets
        with open(path_pathway, 'r') as f1, open(path_conjscan, 'r') as f2:
            pathway_ids = set(line.strip() for line in f1 if line.strip())
            conjscan_ids = set(line.strip() for line in f2 if line.strip())

        overlap = pathway_ids & conjscan_ids

        print(f"Cluster {cluster}:")
        print(f"  Pathway total: {len(pathway_ids)}")
        print(f"  CONJscan total: {len(conjscan_ids)}")
        print(f"  Overlap: {len(overlap)} plasmids")
        print()
        
# Load sequence sizes, removing '>' from header names
sizes_file = os.path.join(output_dir, "sequence_sizes.csv")
sizes_df = pd.read_csv(sizes_file, sep="\t")
sizes_df.columns = ["Plasmid", "Sequence_Length"]
sizes_df["Plasmid"] = sizes_df["Plasmid"].str.lstrip(">")

# Find all pathway cluster files
cluster_files = [f for f in os.listdir(output_dir) if f.startswith("pathway_cluster_") and f.endswith(".txt")]

# Set up plotting
sns.set(style="whitegrid")
for file in sorted(cluster_files, key=lambda x: int(x.split("_")[-1].split(".")[0])):
    cluster_num = file.split("_")[-1].split(".")[0]
    cluster_path = os.path.join(output_dir, file)

    # Load plasmid IDs for this cluster
    with open(cluster_path, 'r') as f:
        plasmid_ids = [line.strip() for line in f if line.strip()]

    # Get matching sizes
    matched = sizes_df[sizes_df["Plasmid"].isin(plasmid_ids)]

    if matched.empty:
        print(f" No matches found for pathway_cluster_{cluster_num}")
        continue

    # Plot
    plt.figure(figsize=(8, 5))
    sns.histplot(matched["Sequence_Length"], bins=30, kde=True, color="skyblue", edgecolor="black")
    plt.title(f"Plasmid Length Distribution - Pathway Cluster {cluster_num}")
    plt.xlabel("Plasmid Length (bp)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.show()

#### Plasmids with integrases

# File path
eggnog_file = r"C:\Users\hayat\Downloads\R_files\data\eggnog_output.emapper.annotations"
output_file = r"C:\Users\hayat\Downloads\R_files\data\plasmids_with_integrases.txt"

plasmids_with_integrase = set()

with open(eggnog_file, 'r') as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split("\t")
        if len(fields) < 21:
            continue
        query, pfam = fields[0], fields[20]
        if "integrase" in pfam.lower() or "rve" in pfam.lower():
            plasmid_id = "_".join(query.split("_")[:-1])
            plasmids_with_integrase.add(plasmid_id)

# Save results
with open(output_file, 'w') as out:
    for plasmid in sorted(plasmids_with_integrase):
        out.write(plasmid + "\n")

print(f"Found {len(plasmids_with_integrase)} plasmids with integrases. Saved to: {output_file}")

# File paths
widespread_file = r"C:\Users\hayat\Downloads\R_files\data\top_abundant_and_widespread_final_plasmid_names.csv"
integrase_file = r"C:\Users\hayat\Downloads\R_files\data\plasmids_with_integrases.txt"
output_file = r"C:\Users\hayat\Downloads\R_files\data\widespread_with_integrase_status.csv"

# Load data
widespread_df = pd.read_csv(widespread_file)
integrase_set = set()

# Load integrase plasmid names
with open(integrase_file, 'r') as f:
    for line in f:
        integrase_set.add(line.strip())

# Assume column is called "Plasmid_ID" or similar — adjust if needed
colname = widespread_df.columns[0]
widespread_df["Integrase"] = widespread_df[colname].apply(lambda x: "+" if x in integrase_set else "-")

# Save to CSV
widespread_df.to_csv(output_file, index=False)
print(f"Output saved to {output_file}")

# Count how many have integrases
n_with_integrase = (widespread_df["Integrase"] == "+").sum()
n_total = len(widespread_df)
print(f"{n_with_integrase} out of {n_total} plasmids have integrases ({n_with_integrase/n_total:.1%})")


status_file = os.path.join(output_dir, "widespread_with_integrase_status.csv")
output_file = os.path.join(output_dir, "widespread_with_integrase_and_length.csv")


# Load widespread plasmid status
status_df = pd.read_csv(status_file)

# Merge on plasmid ID
merged_df = pd.merge(status_df, sizes_df, on="Plasmid", how="left")

# Save result
merged_df.to_csv(output_file, index=False)
print(f" Output saved to: {output_file}")

file_path = r"C:\Users\hayat\Downloads\R_files\data\widespread_with_integrase_and_length.csv"

# Load the data
df = pd.read_csv(file_path)


# Set the style
sns.set(style="whitegrid")

# Create figure
plt.figure(figsize=(10, 6))

# Create a boxplot with jittered stripplot overlay
sns.boxplot(x="Integrase", y="Sequence_Length", data=df, palette="pastel", showfliers=False)
sns.stripplot(x="Integrase", y="Sequence_Length", data=df, jitter=True, color="black", size=3, alpha=0.6)

# Labels and title
plt.xlabel("Integrase Present")
plt.ylabel("Plasmid Length (bp)")
plt.title("Plasmid Length Distribution by Integrase Presence")

# Show the plot
plt.tight_layout()
plt.show()


# Plot percentages with integrases
plt.figure(figsize=(10, 6))
sns.histplot(data=df, x="Sequence_Length", hue="Integrase", bins=30, palette="Set1", kde=False, multiple="stack", alpha=0.7)

# Labels
plt.xlabel("Plasmid Length (bp)")
plt.ylabel("Count")
plt.title("Plasmid Length Distribution Colored by Integrase Presence")
plt.tight_layout()
plt.show()


# Normalize counts to percent within each bin
plt.figure(figsize=(10, 6))
sns.histplot(
    data=df,
    x="Sequence_Length",
    hue="Integrase",
    bins=30,
    palette={"+" : "#d62728", "-" : "#1f77b4"},
    multiple="fill",  # stack and normalize to 100%
    stat="percent"
)

plt.xlabel("Plasmid Length (bp)")
plt.ylabel("Percentage (%)")
plt.title("Relative Plasmid Length Distribution by Integrase Presence")
plt.tight_layout()
plt.show()


