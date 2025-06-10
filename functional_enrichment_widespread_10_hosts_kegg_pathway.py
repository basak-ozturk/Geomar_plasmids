# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 16:43:42 2025

@author: hayat
"""

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
#import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#import re


# Load list of plasmids
with open("C:/Users/hayat/Downloads/R_files/data/widespread_10_more_hosts_list.txt") as f:
    plasmids = set(line.strip() for line in f if line.strip())

# Find the header line (first line starting with '#query')
with open("C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations") as f:
    for i, line in enumerate(f):
        if line.startswith("#query"):
            header_line = i
            break

# Load eggNOG annotations from the correct header line
eggnog_df = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations",
    sep="\t",
    skiprows=header_line,  # skip metadata lines before the actual header
    low_memory=False
)

# Extract plasmid ID from query column
eggnog_df["Plasmid"] = eggnog_df["#query"].str.extract(r"(^[^_]+_[^_]+)")

# Filter only for plasmids in the "15+ hosts" list
filtered = eggnog_df[eggnog_df["Plasmid"].isin(plasmids)]

# Save filtered result
filtered.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/eggnog_widespread_10_hosts.tsv",
    sep="\t",
    index=False
)

def safe_split(pathway):
    if isinstance(pathway, str):
        return [item.strip() for item in pathway.split(',')]
    else:
        return []  
def unify_kegg_id(kegg_id):
    """
    Standardize KEGG pathway IDs to the format 'ko' + 5-digit number with leading zeros.
    Examples:
        'ko1130' -> 'ko01130'
        'map01130' -> 'ko01130'  # if you want to unify prefixes to 'ko'
        'ko00100' -> 'ko00100'
    """
    if not isinstance(kegg_id, str):
        return kegg_id  # return as is if not string

    # Remove known prefixes like 'map', 'ko', etc.
    # If you want to unify all prefixes to 'ko', do so here:
    prefixes = ['ko', 'map', 'rn', 'ec']
    for prefix in prefixes:
        if kegg_id.startswith(prefix):
            num_part = kegg_id[len(prefix):]
            # Pad numeric part to 5 digits
            num_part_padded = num_part.zfill(5)
            return 'ko' + num_part_padded

    # If no known prefix, return as is
    return kegg_id


# First, identify all plasmids in the eggNOG file
eggnog_df["Plasmid"] = eggnog_df["#query"].str.extract(r"(^[^_]+_[^_]+)")

# Create background as plasmids not in foreground
background_df = eggnog_df[~eggnog_df["Plasmid"].isin(plasmids)]

# Remove rows where KEGG_Pathway is "-"
filtered = filtered[filtered["KEGG_Pathway"] != "-"]
background_df = background_df[background_df["KEGG_Pathway"] != "-"]

# Aggregate multiple KEGG pathways per protein
filtered['KEGG_Pathway'] = filtered['KEGG_Pathway'].apply(safe_split)
background_df['KEGG_Pathway'] = background_df['KEGG_Pathway'].apply(safe_split)

# Expand the rows for multiple KEGG pathways
filtered = filtered.explode('KEGG_Pathway')
background_df = background_df.explode('KEGG_Pathway')



# Apply to filtered and background data before counting
filtered['KEGG_Pathway_Unified'] = filtered['KEGG_Pathway'].apply(unify_kegg_id)
background_df['KEGG_Pathway_Unified'] = background_df['KEGG_Pathway'].apply(unify_kegg_id)

# Now count on unified IDs
fg_counts = filtered['KEGG_Pathway_Unified'].value_counts()
bg_counts = background_df['KEGG_Pathway_Unified'].value_counts()

# Proceed with enrichment using fg_counts and bg_counts as before
all_pathways = set(fg_counts.index) | set(bg_counts.index)

results = []
for pathway in all_pathways:
    fg_hits = fg_counts.get(pathway, 0)
    fg_miss = len(filtered) - fg_hits
    bg_hits = bg_counts.get(pathway, 0)
    bg_miss = len(background_df) - bg_hits

    table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    odds_ratio, p = fisher_exact(table, alternative="greater")

    results.append({
        "KEGG_Pathway": pathway,
        "Foreground_Count": fg_hits,
        "Background_Count": bg_hits,
        "Odds_Ratio": odds_ratio,
        "P_Value": p
    })

enrichment_df = pd.DataFrame(results)


# Adjust p-values (optional: Bonferroni or Benjamini-Hochberg)
enrichment_df["P_Value_Adjusted"] = multipletests(enrichment_df["P_Value"], method='fdr_bh')[1]
# Sort by significance
enrichment_df = enrichment_df.sort_values("P_Value")

# Save results
enrichment_df.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/enriched_KEGG_Pathways_10_hosts_vs_all.tsv",
    sep="\t",
    index=False
)

# Set thresholds
significance_threshold = 0.05
odds_ratio_threshold = 1.0  # Only consider enriched terms with OR > 1

# Filter significant terms by adjusted p-value and odds ratio
significant_filtered = enrichment_df[
    (enrichment_df["P_Value_Adjusted"] < significance_threshold) &
    (enrichment_df["Odds_Ratio"] > odds_ratio_threshold)
]

# Summary statistics
summary_stats = {
    'Total_Terms': len(enrichment_df),
    'Significant_Terms_AdjP': len(significant_filtered[significant_filtered["P_Value_Adjusted"] < significance_threshold]),
    'Significant_Terms_AdjP_and_OR': len(significant_filtered)
}

print("Summary statistics:")
for k, v in summary_stats.items():
    print(f"{k}: {v}")

# Print the significant terms and their statistics
print("\nSignificant Enriched Terms:")
for index, row in significant_filtered.iterrows():
    print(f"  KEGG_Pathway: {row['KEGG_Pathway']}")
    print(f"    Adjusted P-value: {row['P_Value_Adjusted']:.3e}")
    print(f"    Odds Ratio: {row['Odds_Ratio']:.3f}")
    print(f"    Foreground Count: {row['Foreground_Count']}")
    print(f"    Background Count: {row['Background_Count']}")
    print("-" * 40)
    
# --- Visualization ---
# Bar Plot of Odds Ratios

# Set up the figure and axes
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

# Sort the DataFrame by Odds Ratio for better visualization
significant_filtered_sorted = significant_filtered.sort_values(
    "Odds_Ratio", ascending=False
)

# Create the bar plot
barplot = sns.barplot(
    x="Odds_Ratio",
    y="KEGG_Pathway",
    data=significant_filtered_sorted,
    palette="viridis",
)

# Add labels and title
plt.xlabel("Odds Ratio", fontsize=12)
plt.ylabel("KEGG Pathway", fontsize=12)
plt.title("Odds Ratio of Significant KEGG Pathways", fontsize=14)

# # Add annotations (values) to the bars
# for p in barplot.patches:
#     width = p.get_width()
#     plt.text(
#         5,
#         p.get_y() + p.get_height() / 2,
#         f"{width:.1f}",
#         ha="left",
#         va="center",
#         fontsize=10,
#         color="black",
#     )

# Show the plot
plt.tight_layout()
plt.show()



