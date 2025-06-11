# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 10:45:33 2025

@author: hayat
"""

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# Load list of plasmids
with open("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid_names.csv") as f:
    plasmids = set(line.strip() for line in f if line.strip())

# Find header line in eggNOG annotations
with open("C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations") as f:
    for i, line in enumerate(f):
        if line.startswith("#query"):
            header_line = i
            break

# Load eggNOG annotations
eggnog_df = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations",
    sep="\t",
    skiprows=header_line,
    low_memory=False
)

# Extract plasmid IDs once
eggnog_df["Plasmid"] = eggnog_df["#query"].str.extract(r"(^[^_]+_[^_]+)")

# Filter foreground and background
foreground = eggnog_df[eggnog_df["Plasmid"].isin(plasmids)].copy()
background = eggnog_df[~eggnog_df["Plasmid"].isin(plasmids)].copy()

# Remove rows with missing PFAMs
foreground = foreground[foreground["PFAMs"] != "-"]
background = background[background["PFAMs"] != "-"]

# Split multi-PFAM annotations into lists
foreground['PFAMs_list'] = foreground['PFAMs'].str.split(',')
background['PFAMs_list'] = background['PFAMs'].str.split(',')

# Explode lists so each PFAM has its own row
foreground_exp = foreground.explode('PFAMs_list')
background_exp = background.explode('PFAMs_list')

# Strip whitespace
foreground_exp['PFAMs_list'] = foreground_exp['PFAMs_list'].str.strip()
background_exp['PFAMs_list'] = background_exp['PFAMs_list'].str.strip()

# Define PFAM groups to merge related PFAMs into categories
pfam_group_map = {
    'Arm-DNA-bind_4': 'Phage integrases',
    'Phage_int_SAM_3': 'Phage integrases',
    'Phage_integrase': 'Phage integrases',
    'Arm-DNA-bind_3': 'Phage integrases',
    'GFO_IDH_MocA': 'GFO_IDH_MocA_C',
    'GFO_IDH_MocA_C': 'GFO_IDH_MocA_C',
    # Add other groups here as needed
}

def map_to_group(pfam):
    return pfam_group_map.get(pfam, pfam)  # default to original PFAM if no group

# Map PFAMs to groups
foreground_exp['PFAM_Group'] = foreground_exp['PFAMs_list'].apply(map_to_group)
background_exp['PFAM_Group'] = background_exp['PFAMs_list'].apply(map_to_group)

# Count occurrences by PFAM group
fg_counts = foreground_exp['PFAM_Group'].value_counts()
bg_counts = background_exp['PFAM_Group'].value_counts()

all_pfams = fg_counts.index.union(bg_counts.index)

fg_total = len(foreground_exp)
bg_total = len(background_exp)

# Enrichment analysis
results = []
for pfam in all_pfams:
    fg_hits = fg_counts.get(pfam, 0)
    fg_miss = fg_total - fg_hits
    bg_hits = bg_counts.get(pfam, 0)
    bg_miss = bg_total - bg_hits

    contingency_table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    odds_ratio, p_value = fisher_exact(contingency_table, alternative="greater")

    results.append({
        "PFAMs": pfam,
        "Foreground_Count": fg_hits,
        "Background_Count": bg_hits,
        "Odds_Ratio": odds_ratio,
        "P_Value": p_value
    })

enrichment_df = pd.DataFrame(results)

# Multiple testing correction
enrichment_df["P_Value_Adjusted"] = multipletests(enrichment_df["P_Value"], method='fdr_bh')[1]

# Set thresholds
significance_threshold = 0.05
odds_ratio_threshold = 1.0

# Dynamic minimum foreground count cutoff (e.g., 10th percentile)
min_fg_count = 10
print(f"Using minimum foreground count cutoff: {min_fg_count}")

# Filter significant terms with cutoff
significant_filtered = enrichment_df[
    (enrichment_df["P_Value_Adjusted"] < significance_threshold) &
    (enrichment_df["Odds_Ratio"] > odds_ratio_threshold) &
    (enrichment_df["Foreground_Count"] >= min_fg_count)
].copy()

# Sort for plotting
significant_filtered.sort_values("Odds_Ratio", ascending=False, inplace=True)

# Summary
print("Summary statistics:")
print(f"Total PFAM groups tested: {len(enrichment_df)}")
print(f"Significant PFAM groups (Adj P < {significance_threshold}): {len(enrichment_df[enrichment_df['P_Value_Adjusted'] < significance_threshold])}")
print(f"Significant PFAM groups with OR > {odds_ratio_threshold} and FG count >= {min_fg_count}: {len(significant_filtered)}")

print("\nSignificant Enriched PFAM Domains (grouped):")
for _, row in significant_filtered.iterrows():
    print(f"{row['PFAMs']}: Adj P={row['P_Value_Adjusted']:.3e}, OR={row['Odds_Ratio']:.2f}, FG={row['Foreground_Count']}, BG={row['Background_Count']}")

# Visualization
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

barplot = sns.barplot(
    x="Odds_Ratio",
    y="PFAMs",
    data=significant_filtered,
    palette="viridis"
)

plt.xlabel("Odds Ratio", fontsize=12)
plt.ylabel("PFAM Domain Group", fontsize=12)
plt.title("Odds Ratio of Significant PFAM Domain Groups", fontsize=14)

# Dynamic annotation positioning
for p in barplot.patches:
    width = p.get_width()
    barplot.text(
        width + 0.1,
        p.get_y() + p.get_height() / 2,
        f"{width:.2f}",
        ha="left",
        va="center",
        fontsize=10,
        color="black"
    )

plt.tight_layout()
plt.show()

# Save filtered significant enrichment results
significant_filtered.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/enriched_PFAMs_all_widespread_significant.tsv",
    sep="\t",
    index=False
)