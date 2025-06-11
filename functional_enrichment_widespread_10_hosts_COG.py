# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 10:21:29 2025

@author: hayat
"""

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns


with open("C:/Users/hayat/Downloads/R_files/data/widespread_10_more_hosts_list.txt") as f:
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



# Filter out unknown or uninformative categories first (including NaNs)
exclude_categories = {'S', 'R', 'X', '-', None, ''}
filtered_df = eggnog_df[~eggnog_df['COG_category'].isin(exclude_categories)].copy()

# Extract plasmid IDs
filtered_df["Plasmid"] = filtered_df["#query"].str.extract(r"(^[^_]+_[^_]+)")

# Safely split concatenated COG category strings into list of single letters
def split_cog_categories(cog):
    if isinstance(cog, str):
        return list(cog.strip())
    else:
        return []

filtered_df['COG_category'] = filtered_df['COG_category'].apply(split_cog_categories)

# Explode into separate rows per COG category letter
filtered_df = filtered_df.explode('COG_category')

# Strip whitespace just in case
filtered_df['COG_category'] = filtered_df['COG_category'].str.strip()

# Now define foreground and background
foreground = filtered_df[filtered_df["Plasmid"].isin(plasmids)].copy()
background = filtered_df[~filtered_df["Plasmid"].isin(plasmids)].copy()

# Count occurrences of each COG category
fg_counts = foreground['COG_category'].value_counts()
bg_counts = background['COG_category'].value_counts()

all_cogs = fg_counts.index.union(bg_counts.index)

fg_total = len(foreground)
bg_total = len(background)


# Proceed with enrichment as before
results = []
for cog in all_cogs:
    fg_hits = fg_counts.get(cog, 0)
    fg_miss = fg_total - fg_hits
    bg_hits = bg_counts.get(cog, 0)
    bg_miss = bg_total - bg_hits

    table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    odds_ratio, p_value = fisher_exact(table, alternative='greater')

    results.append({
        'COG_category': cog,
        'Foreground_Count': fg_hits,
        'Background_Count': bg_hits,
        'Odds_Ratio': odds_ratio,
        'P_Value': p_value
    })

enrichment_df = pd.DataFrame(results)


# Adjust p-values for multiple testing
enrichment_df['P_Value_Adjusted'] = multipletests(enrichment_df['P_Value'], method='fdr_bh')[1]

# Optional: dynamic minimum foreground count cutoff (e.g., 10th percentile)
min_fg_count = int(enrichment_df['Foreground_Count'].quantile(0.10))
print(f"Using minimum foreground count cutoff: {min_fg_count}")

# Filter significant enriched COG categories
significant_cogs = enrichment_df[
    (enrichment_df['P_Value_Adjusted'] < 0.05) &
    (enrichment_df['Odds_Ratio'] > 1) &
    (enrichment_df['Foreground_Count'] >= min_fg_count)
].sort_values('Odds_Ratio', ascending=False)

print("Significant enriched COG categories:")
print(significant_cogs)

plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")

barplot = sns.barplot(
    x="Odds_Ratio",
    y="COG_category",
    data=significant_cogs,
    palette="viridis"
)

plt.xlabel("Odds Ratio", fontsize=14)
plt.ylabel("COG category", fontsize=14)
plt.title("Significant Enriched COG categories", fontsize=16)

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
