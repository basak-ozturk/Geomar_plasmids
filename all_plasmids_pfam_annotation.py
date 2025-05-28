# -*- coding: utf-8 -*-
"""
Created on Wed May 28 15:35:34 2025

@author: hayat
"""

import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

input_file = "C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations"
pfam_counter = Counter()

with open(input_file, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) < 1:
            continue
        pfam_col = fields[-1].strip()
        if pfam_col and pfam_col != "-":
            pfams = [pfam.strip() for pfam in pfam_col.split(",") if pfam.strip()]
            pfam_counter.update(pfams)

# Convert to DataFrame
top_pfam_df = pd.DataFrame(pfam_counter.most_common(100), columns=["PFAM", "Count"])

# Save to CSV
output_file = "C:/Users/hayat/Downloads/R_files/data/top_100_pfams_in_all_annotations.csv"
top_pfam_df.to_csv(output_file, index=False)

print(f"Top 100 PFAMs written to: {output_file}")

top_pfam_df_curated = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_pfam_all_plasmids_for_curation_curated_consolidated.csv")

# Optional: Sort by Count for consistent bar order
top_pfam_df_curated = top_pfam_df_curated.sort_values("Count", ascending=False)

# Assign colors by Category
plt.figure(figsize=(14, 10))
sns.set_style("whitegrid")

# Use Seaborn barplot to color by Category
ax = sns.barplot(
    data=top_pfam_df_curated,
    y="PFAM", x="Count",
    hue="Category",
    dodge=False,
    palette="Paired"  
)

# Adjust aesthetics
plt.title("Top PFAM Domains in All Plasmid Proteins (Categorized)", fontsize=16)
plt.xlabel("Count", fontsize=12)
plt.ylabel("PFAM", fontsize=12)
plt.yticks(fontsize=8)
plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/all_plasmids_top_pfam_plot_categorized.png", dpi=300)
plt.show()