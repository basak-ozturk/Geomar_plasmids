# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 14:00:56 2025

@author: hayat
"""

from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Define input files
file1 = "C:/Users/hayat/Downloads/R_files/data/plasmids_with_plasmid_hallmark.txt"
file2 = "C:/Users/hayat/Downloads/R_files/data/oriT_alloriT_blast_results_names.txt"
file3 = "C:/Users/hayat/Downloads/R_files/data/unique_CONJscan_hit_ids.txt"

# Read content into sets
with open(file1, "r") as f:
    set1 = set(line.strip() for line in f if line.strip())

with open(file2, "r") as f:
    set2 = set(line.strip() for line in f if line.strip())

with open(file3, "r") as f:
    set3 = set(line.strip() for line in f if line.strip())

# Create the Venn diagram
plt.figure(figsize=(8, 8))
venn = venn3(
    [set1, set2, set3],
    set_labels=("Plasmid Hallmarks", r"$\it{oriT}$", "CONJscan"),
    set_colors=("#E63946", "#457B9D", "#2A9D8F"),  # nicer colors
    alpha=0.6
)

plt.title("Venn Diagram of PCUs with Plasmid Hallmarks, CONJscan hits, and r$\it{oriT}$", fontsize=16, fontweight='bold')

# Save as SVG
output_path = "C:/Users/hayat/Downloads/R_files/graphs/venn_diagram_orit_conjscan_hallmark.svg"
plt.savefig(output_path, format="svg")
plt.show()
plt.close()

print(f"Venn diagram saved to: {output_path}")
