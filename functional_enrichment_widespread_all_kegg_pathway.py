# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 16:43:42 2025

@author: hayat
"""

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# --- Helper functions ---

def safe_split(pathway):
    """Split KEGG pathway string by commas safely, handle missing values."""
    if isinstance(pathway, str):
        return [item.strip() for item in pathway.split(',')]
    return []

def unify_kegg_id(kegg_id):
    """
    Standardize KEGG pathway IDs to format 'ko' + 5-digit number with leading zeros.
    Unify prefixes like 'map', 'rn', 'ec' to 'ko'.
    """
    if not isinstance(kegg_id, str):
        return kegg_id
    prefixes = ['ko', 'map', 'rn', 'ec']
    for prefix in prefixes:
        if kegg_id.startswith(prefix):
            num_part = kegg_id[len(prefix):]
            num_part_padded = num_part.zfill(5)
            return 'ko' + num_part_padded
    return kegg_id

# --- Load data ---

# Load plasmid list
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

# Filter foreground and background plasmids
foreground = eggnog_df[eggnog_df["Plasmid"].isin(plasmids)].copy()
background = eggnog_df[~eggnog_df["Plasmid"].isin(plasmids)].copy()

# Remove rows with missing KEGG_Pathway
foreground = foreground[foreground["KEGG_Pathway"] != "-"]
background = background[background["KEGG_Pathway"] != "-"]

# Split multiple KEGG pathways into lists
foreground["KEGG_Pathway"] = foreground["KEGG_Pathway"].apply(safe_split)
background["KEGG_Pathway"] = background["KEGG_Pathway"].apply(safe_split)

# Explode lists into separate rows
foreground = foreground.explode("KEGG_Pathway")
background = background.explode("KEGG_Pathway")

# Standardize KEGG pathway IDs
foreground["KEGG_Pathway_Unified"] = foreground["KEGG_Pathway"].apply(unify_kegg_id)
background["KEGG_Pathway_Unified"] = background["KEGG_Pathway"].apply(unify_kegg_id)

# Count occurrences per pathway
fg_counts = foreground["KEGG_Pathway_Unified"].value_counts()
bg_counts = background["KEGG_Pathway_Unified"].value_counts()

# All unique pathways
all_pathways = fg_counts.index.union(bg_counts.index)

# --- Enrichment analysis ---

results = []
fg_total = len(foreground)
bg_total = len(background)

for pathway in all_pathways:
    fg_hits = fg_counts.get(pathway, 0)
    fg_miss = fg_total - fg_hits
    bg_hits = bg_counts.get(pathway, 0)
    bg_miss = bg_total - bg_hits

    contingency_table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    odds_ratio, p_value = fisher_exact(contingency_table, alternative="greater")

    results.append({
        "KEGG_Pathway": pathway,
        "Foreground_Count": fg_hits,
        "Background_Count": bg_hits,
        "Odds_Ratio": odds_ratio,
        "P_Value": p_value
    })

enrichment_df = pd.DataFrame(results)

# Adjust p-values for multiple testing (FDR)
enrichment_df["P_Value_Adjusted"] = multipletests(enrichment_df["P_Value"], method="fdr_bh")[1]

# Filter significant terms with thresholds
significance_threshold = 0.05
odds_ratio_threshold = 1.0
min_fg_count = 10

significant_filtered = enrichment_df[
    (enrichment_df["P_Value_Adjusted"] < significance_threshold) &
    (enrichment_df["Odds_Ratio"] > odds_ratio_threshold) &
    (enrichment_df["Foreground_Count"] >= min_fg_count)
].copy()

# Sort by Odds Ratio descending for visualization
significant_filtered.sort_values("Odds_Ratio", ascending=False, inplace=True)

# --- Output summary ---

print("Summary statistics:")
print(f"Total terms tested: {len(enrichment_df)}")
print(f"Significant terms (Adj P < {significance_threshold}): {len(enrichment_df[enrichment_df['P_Value_Adjusted'] < significance_threshold])}")
print(f"Significant terms with OR > {odds_ratio_threshold} and FG count >= {min_fg_count}: {len(significant_filtered)}")

print("\nSignificant enriched KEGG pathways:")
for _, row in significant_filtered.iterrows():
    print(f"{row['KEGG_Pathway']}: Adj P={row['P_Value_Adjusted']:.3e}, OR={row['Odds_Ratio']:.2f}, FG={row['Foreground_Count']}, BG={row['Background_Count']}")

# --- Visualization ---

plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")

barplot = sns.barplot(
    x="Odds_Ratio",
    y="KEGG_Pathway",
    data=significant_filtered,
    palette="viridis"
)

plt.xlabel("Odds Ratio", fontsize=14)
plt.ylabel("KEGG Pathway", fontsize=14)
plt.title("Significant Enriched KEGG Pathways", fontsize=16)

# Annotate bars with odds ratio values
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
    "C:/Users/hayat/Downloads/R_files/data/enriched_KEGG_Pathways_all_widespread_significant.tsv",
    sep="\t",
    index=False
)
