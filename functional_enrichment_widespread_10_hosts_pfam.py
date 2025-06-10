import pandas as pd

from scipy.stats import fisher_exact

from statsmodels.stats.multitest import multipletests

import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

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


# First, identify all plasmids in the eggNOG file
eggnog_df["Plasmid"] = eggnog_df["#query"].str.extract(r"(^[^_]+_[^_]+)")

# Create background as plasmids not in foreground
background_df = eggnog_df[~eggnog_df["Plasmid"].isin(plasmids)]

# Remove rows where PFAMs is "-"
filtered = filtered[filtered["PFAMs"] != "-"]
background_df = background_df[background_df["PFAMs"] != "-"]

# Now compute counts
fg_counts = filtered['PFAMs'].value_counts()
bg_counts = background_df['PFAMs'].value_counts()

fg_total = len(filtered)
bg_total = len(background_df)



# Store enrichment results # Union of all PFAMss (update this after defining fg_counts and bg_counts)
all_PFAMss = set(fg_counts.index) | set(bg_counts.index)

results = []

for desc in all_PFAMss:
    fg_hits = fg_counts.get(desc, 0)
    fg_miss = fg_total - fg_hits
    bg_hits = bg_counts.get(desc, 0)
    bg_miss = bg_total - bg_hits

    # 2x2 contingency table
    table = [[fg_hits, fg_miss], [bg_hits, bg_miss]]

    # Fisher’s exact test (right-sided = overrepresentation in foreground)
    odds_ratio, p = fisher_exact(table, alternative="greater")

    results.append({
        "PFAMs": desc,
        "Foreground_Count": fg_hits,
        "Background_Count": bg_hits,
        "Odds_Ratio": odds_ratio,
        "P_Value": p
    })

# Turn into dataframe
enrichment_df = pd.DataFrame(results)

# Adjust p-values (optional: Bonferroni or Benjamini-Hochberg)
enrichment_df["P_Value_Adjusted"] = multipletests(enrichment_df["P_Value"], method='fdr_bh')[1]
# Sort by significance
enrichment_df = enrichment_df.sort_values("P_Value")

# Save results
enrichment_df.to_csv(
    "C:/Users/hayat/Downloads/R_files/data/enriched_PFAMss_10_hosts_vs_all.tsv",
    sep="\t",
    index=False
)

enrichment_df['log2_OR'] = np.log2(enrichment_df['Odds_Ratio'].replace(0, np.nan))
enrichment_df['neg_log10_p'] = -np.log10(enrichment_df['P_Value_Adjusted'])

# After creating enrichment_df and adding P_Value_Adjusted:
enrichment_df["P_Value_Adjusted"] = multipletests(enrichment_df["P_Value"], method='fdr_bh')[1]

# Set thresholds
significance_threshold = 0.05
odds_ratio_threshold = 1.0  # Only consider enriched terms with OR > 1

# Filter significant terms by adjusted p-value
significant_terms = enrichment_df[enrichment_df["P_Value_Adjusted"] < significance_threshold]

# Further filter by odds ratio threshold
significant_filtered = significant_terms[significant_terms["Odds_Ratio"] > odds_ratio_threshold]

# Summary statistics
summary_stats = {
    'Total_Terms': len(enrichment_df),
    'Significant_Terms_AdjP': len(significant_terms),
    'Significant_Terms_AdjP_and_OR': len(significant_filtered)
}

print("Summary statistics:")
for k, v in summary_stats.items():
    print(f"{k}: {v}")

# Relaxed thresholds for exploratory purposes
exploratory_significance_threshold = 0.1
exploratory_odds_ratio_threshold = 1.0

exploratory_filtered = enrichment_df[
    (enrichment_df["P_Value_Adjusted"] < exploratory_significance_threshold) &
    (enrichment_df["Odds_Ratio"] > exploratory_odds_ratio_threshold)
]

print(f"Terms passing relaxed thresholds: {len(exploratory_filtered)}")

# Filter for significant terms
significant_filtered = enrichment_df[
    (enrichment_df["P_Value_Adjusted"] < significance_threshold) &
    (enrichment_df["Odds_Ratio"] > odds_ratio_threshold)
]

# Print the significant terms and their statistics
print("\nSignificant Enriched Terms:")
for index, row in significant_filtered.iterrows():
    print(f"  PFAMs: {row['PFAMs']}")
    print(f"    Adjusted P-value: {row['P_Value_Adjusted']:.3e}")
    print(f"    Odds Ratio: {row['Odds_Ratio']:.3f}")
    print(f"    Foreground Count: {row['Foreground_Count']}")
    print(f"    Background Count: {row['Background_Count']}")
    print("-" * 40)
    
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
    y="PFAMs",
    data=significant_filtered_sorted,
    palette="viridis",
)

# Add labels and title
plt.xlabel("Odds Ratio", fontsize=12)
plt.ylabel("PFAM Domain", fontsize=12)
plt.title("Odds Ratio of Significant PFAM Domains", fontsize=14)

# Add annotations (values) to the bars
for p in barplot.patches:
    width = p.get_width()
    plt.text(
        5,
        p.get_y() + p.get_height() / 2,
        f"{width:.1f}",
        ha="left",
        va="center",
        fontsize=10,
        color="black",
    )

# Show the plot
plt.tight_layout()
plt.show()
