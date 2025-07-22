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
import numpy as np

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
Geodia_plasmid_file = "C:/Users/hayat/Downloads/R_files/data/plasmids_in_Agelas.tsv"
Geodia_df = pd.read_csv(Geodia_plasmid_file, sep="\t")
plasmids = set(Geodia_df.iloc[:, 0])  # First column assumed to be plasmid names

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
min_fg_count = 20

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


# --- Add log2(OR) and 95% CI to the DataFrame ---

def compute_log2_ci(row):
    a = row["Foreground_Count"]
    b = fg_total - a
    c = row["Background_Count"]
    d = bg_total - c

    # Avoid zero counts for stability
    a = max(a, 1)
    b = max(b, 1)
    c = max(c, 1)
    d = max(d, 1)

    log_or = np.log2(row["Odds_Ratio"])
    se = np.sqrt(1/a + 1/b + 1/c + 1/d) / np.log(2)  # convert SE to log2 scale
    ci_lower = log_or - 1.96 * se
    ci_upper = log_or + 1.96 * se
    return pd.Series([log_or, ci_lower, ci_upper])

significant_filtered[["log2_OR", "log2_CI_Lower", "log2_CI_Upper"]] = significant_filtered.apply(compute_log2_ci, axis=1)

# --- Forest plot ---

# Load annotation CSV
annotations = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/agelas_enriched_pathways_annotations.csv", header=None, names=["ID", "Label"])

# Remove whitespace
annotations["ID"] = annotations["ID"].str.strip()
annotations["Label"] = annotations["Label"].str.strip()

# Merge with the main dataframe
significant_filtered = significant_filtered.merge(
    annotations,
    left_on="KEGG_Pathway",
    right_on="ID",
    how="left"
)

significant_filtered["Y_Label"] = significant_filtered.apply(
    lambda row: f"{row['Label']} ({row['KEGG_Pathway']})" if pd.notna(row["Label"]) else row["KEGG_Pathway"],
    axis=1
)

# Sort for better layout
significant_filtered = significant_filtered.sort_values("log2_OR", ascending=True)

plt.figure(figsize=(10, 8))
sns.set(style="whitegrid")

plt.errorbar(
    x=significant_filtered["log2_OR"],
    y=significant_filtered["Y_Label"],
    xerr=[
        significant_filtered["log2_OR"] - significant_filtered["log2_CI_Lower"],
        significant_filtered["log2_CI_Upper"] - significant_filtered["log2_OR"]
    ],
    fmt='o',
    color='black',
    ecolor='gray',
    elinewidth=2,
    capsize=4
)

plt.axvline(0, color='red', linestyle='--', label='OR = 1')
plt.xlabel(r'$\log_2(\mathrm{Odds\ Ratio})$', fontsize=12)
plt.ylabel("Pathway", fontsize=12)
plt.title(r'Pathway Enrichment in Agelas Plasmids' + '\n' +
          r'$\left(\log_2(\mathrm{Odds\ Ratio\ with\ 95\%\ CI})\right)$')

plt.tight_layout()
plt.legend()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/agelas_enriched_pathways_OR_with_CI.png", dpi=300)
plt.show()


# Save filtered significant enrichment results
significant_filtered.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/enriched_KEGG_Pathways_Agelas.tsv",
    sep="\t",
    index=False
)
