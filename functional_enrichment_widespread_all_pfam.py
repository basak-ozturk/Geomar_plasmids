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
from statsmodels.stats.contingency_tables import Table2x2
import numpy as np
# Load list of plasmids
with open("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_names.txt") as f:
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

results = []
for pfam in all_pfams:
    fg_hits = fg_counts.get(pfam, 0)
    fg_miss = fg_total - fg_hits
    bg_hits = bg_counts.get(pfam, 0)
    bg_miss = bg_total - bg_hits

    contingency_table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    table = Table2x2(contingency_table)

    odds_ratio = table.oddsratio
    ci_low, ci_upp = table.oddsratio_confint()
    p_value = fisher_exact(contingency_table, alternative="greater")[1]

    results.append({
        "PFAMs": pfam,
        "Foreground_Count": fg_hits,
        "Background_Count": bg_hits,
        "Odds_Ratio": odds_ratio,
        "CI_Lower": ci_low,
        "CI_Upper": ci_upp,
        "P_Value": p_value
    })

enrichment_df = pd.DataFrame(results)
enrichment_df["P_Value_Adjusted"] = multipletests(enrichment_df["P_Value"], method='fdr_bh')[1]


# Set thresholds
significance_threshold = 0.05
odds_ratio_threshold = 1.0


min_fg_count = 5
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

# Filter again
significant_filtered = enrichment_df[
    (enrichment_df["P_Value_Adjusted"] < significance_threshold) &
    (enrichment_df["Odds_Ratio"] > odds_ratio_threshold) &
    (enrichment_df["Foreground_Count"] >= min_fg_count)
].copy()

# Add log2-transformed odds ratio and CI bounds
significant_filtered["log2_OR"] = np.log2(significant_filtered["Odds_Ratio"])
significant_filtered["log2_CI_Lower"] = np.log2(significant_filtered["CI_Lower"])
significant_filtered["log2_CI_Upper"] = np.log2(significant_filtered["CI_Upper"])

# Sort for plotting
significant_filtered.sort_values("log2_OR", ascending=True, inplace=True)

# Plot with error bars (forest plot style)
plt.figure(figsize=(10, 8))
sns.set(style="whitegrid")

plt.errorbar(
    x=significant_filtered["log2_OR"],
    y=significant_filtered["PFAMs"],
    xerr=[significant_filtered["log2_OR"] - significant_filtered["log2_CI_Lower"],
          significant_filtered["log2_CI_Upper"] - significant_filtered["log2_OR"]],
    fmt='o',
    color='black',
    ecolor='gray',
    elinewidth=2,
    capsize=4
)

plt.axvline(0, color='red', linestyle='--', label='OR = 1')
plt.xlabel(r'$\log_2(\mathrm{Odds\ Ratio})$', fontsize=12)
plt.ylabel("PFAM Domain Group", fontsize=12)
plt.title(r'PFAM Enrichment in Abundant and Widespread Plasmids' + '\n' +
          r'$\left(\log_2(\mathrm{Odds\ Ratio\ with\ 95\%\ CI})\right)$')

plt.tight_layout()
plt.legend()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/enriched_PFAMs_OR_with_CI.png", dpi=300)

plt.show()

# Save filtered significant enrichment results
significant_filtered.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/enriched_PFAMs_all_widespread_significant_with_CI.tsv",
    sep="\t",
    index=False
)

#foreground.to_csv("C:/Users/hayat/Downloads/R_files/data/filtered_foreground_eggnog_annotations.tsv", sep="\t", index=False)

