import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
#import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load list of plasmids
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

# Extract plasmid IDs once
eggnog_df["Plasmid"] = eggnog_df["#query"].str.extract(r"(^[^_]+_[^_]+)")

# Filter foreground and background
foreground = eggnog_df[eggnog_df["Plasmid"].isin(plasmids)].copy()
background = eggnog_df[~eggnog_df["Plasmid"].isin(plasmids)].copy()

# Remove rows with missing PFAMs
foreground = foreground[foreground["PFAMs"] != "-"]
background = background[background["PFAMs"] != "-"]

# Count PFAM occurrences
fg_counts = foreground["PFAMs"].value_counts()
bg_counts = background["PFAMs"].value_counts()

# All unique PFAMs
all_pfams = fg_counts.index.union(bg_counts.index)

# Total counts
fg_total = len(foreground)
bg_total = len(background)

# Enrichment analysis
results = []
for pfam in all_pfams:
    fg_hits = fg_counts.get(pfam, 0)
    fg_miss = fg_total - fg_hits
    bg_hits = bg_counts.get(pfam, 0)
    bg_miss = bg_total - bg_hits

    table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]
    odds_ratio, p_value = fisher_exact(table, alternative="greater")

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

# Thresholds
significance_threshold = 0.05
odds_ratio_threshold = 1.0
min_fg_count = 8  # Minimum foreground count cutoff

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
print(f"Total PFAMs tested: {len(enrichment_df)}")
print(f"Significant PFAMs (Adj P < {significance_threshold}): {len(enrichment_df[enrichment_df['P_Value_Adjusted'] < significance_threshold])}")
print(f"Significant PFAMs with OR > {odds_ratio_threshold} and FG count >= {min_fg_count}: {len(significant_filtered)}")

print("\nSignificant Enriched PFAM Domains:")
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
plt.ylabel("PFAM Domain", fontsize=12)
plt.title("Odds Ratio of Significant PFAM Domains", fontsize=14)

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
