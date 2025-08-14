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
cog_counter = Counter()

with open(input_file, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split('\t')

        if len(fields) < 1:
            continue

        # Adjust column indices if necessary based on your file format
        pfam_col = fields[-1].strip()  # last column for PFAM
        cog_col = fields[6].strip()   

        # PFAM extraction
        if pfam_col and pfam_col != "-":
            pfams = [pfam.strip() for pfam in pfam_col.split(",") if pfam.strip()]
            pfam_counter.update(pfams)

        # COG extraction
        if cog_col and cog_col != "-":
            cogs = list(cog_col)
            cogs = [c for c in cogs if c not in ["-", "S", "R", ]]
            cog_counter.update(cogs)

# PFAM counts → DataFrame
top_pfam_df = pd.DataFrame(pfam_counter.most_common(100), columns=["PFAM", "Count"])
#top_pfam_df.to_csv("C:/Users/hayat/Downloads/R_files/data/top_50_pfams_in_all_annotations.csv", index=False)

print("Top 50 PFAMs written to CSV.")

# Plot curated PFAM
top_pfam_df_curated = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_pfam_all_plasmids_for_curation_curated_consolidated.csv")
top_pfam_df_curated = top_pfam_df_curated.sort_values("Count", ascending=False)

category_palette = {
    "Toxin-antitoxin": "#D55E00",               # reddish-orange
    "Oxidoreductase": "#0072B2",                # blue
    "Signalling": "#009E73",                    # bluish green
    "Broad/uncharacterized function": "#999999",# gray
    "Restriction-DNA modification": "#E69F00",  # orange
    "Conjugation": "#56B4E9",                   # light blue
    "Transposase": "#F0E442",                   # yellow
    "Integrase": "#CC79A7",                     # magenta
    "Regulatory protein": "#882255",            # deep burgundy
    "Transporter": "#44AA99",                   # teal
    "Secretion system": "#117733",               # dark green
    "Peptidase": "#AA4499"  # deep violet / plum

}



plt.figure(figsize=(14, 10))
sns.set_style("whitegrid")

ax = sns.barplot(
    data=top_pfam_df_curated,
    y="PFAM", x="Count",
    hue="Category",
    dodge=False,
    palette=category_palette

)

plt.title("Top PFAM Domains in All Plasmid Proteins (Categorized)", fontsize=16)
plt.xlabel("Count", fontsize=12)
plt.ylabel("PFAM", fontsize=12)
plt.yticks(fontsize=8)
plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/all_plasmids_top_pfam_plot_categorized.png", dpi=300)
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/all_plasmids_top_pfam_plot_categorized.svg", dpi=300)

plt.show()

# COG counts → DataFrame
cog_counts = pd.DataFrame(cog_counter.items(), columns=["COG", "Count"])
cog_counts = cog_counts.sort_values("Count", ascending=False)

# Map COG to functional names
cog_function_map = {
    "J": "Translation (J)",
    "A": "RNA processing and modification (A)",
    "K": "Transcription (K)",
    "L": "Replication, recombination and repair (L)",
    "B": "Chromatin structure and dynamics (B)",
    "D": "Cell cycle control (D)",
    "Y": "Nuclear structure (Y)",
    "V": "Defense mechanisms (V)",
    "T": "Signal transduction (T)",
    "M": "Cell wall/membrane biogenesis (M)",
    "N": "Cell motility (N)",
    "U": "Intracellular trafficking (U)",
    "O": "Posttranslational modification (O)",
    "C": "Energy production (C)",
    "G": "Carbohydrate metabolism (G)",
    "E": "Amino acid metabolism (E)",
    "F": "Nucleotide metabolism (F)",
    "H": "Coenzyme metabolism (H)",
    "I": "Lipid metabolism (I)",
    "P": "Inorganic ion transport (P)",
    "Q": "Secondary metabolites biosynthesis (Q)",
    "R": "General function prediction (R)",
}

cog_counts["Function"] = cog_counts["COG"].map(cog_function_map)

# COG barplot
plt.figure(figsize=(10, 8))
sns.barplot(data=cog_counts, y="Function", x="Count", color="skyblue")
plt.xlabel("Number of Proteins")
plt.ylabel("COG Functional Category")
plt.title("COG Categories in All Plasmid Proteins")
plt.tight_layout()
#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/all_plasmids_cog_category_plot_filtered.png", dpi=300)
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/all_plasmids_cog_category_plot_filtered.svg", dpi=300)

plt.show()
