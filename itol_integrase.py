# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 09:02:08 2025

@author: hayat
"""

with open("C:/Users/hayat/Downloads/R_files/data/integrase_proteins_widepsread.txt") as f:
    Integrase_ids = [line.strip().rsplit('_', 1)[0] for line in f]

Integrase_plasmids = set(Integrase_ids)


# Read all plasmid names
with open(r"C:\Users\hayat\Downloads\R_files\data\top_abundant_and_widespread_plasmid_names.txt") as f:
    all_plasmids = [line.strip() for line in f if line.strip()]

# Create iTOL colorstrip dataset content
output_lines = [
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tIntegrase",
    "COLOR\t#000000",
    "LEGEND_TITLE\tIntegrase",
    "LEGEND_SHAPES\t1\t1",
    "LEGEND_COLORS\t#f77f00\t#274e13",
    "LEGEND_LABELS\tIntegrase+\tIntegrase-",
    "DATA"
]

# Assign colors and labels
for plasmid in sorted(all_plasmids):
    if plasmid in Integrase_plasmids:
        color = "#f77f00"
        label = "Integrase+"
    else:
        color = "#274e13"
        label = "Integrase-"
    output_lines.append(f"{plasmid}\t{color}\t{label}")

# Save to file
with open(r"C:\Users\hayat\Downloads\R_files\data\itol_integrase.txt", "w") as out:
    out.write("\n".join(output_lines))

print("iTOL file saved as 'itol_Integrase.txt'")
