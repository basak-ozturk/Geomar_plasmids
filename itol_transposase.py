# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 09:02:08 2025

@author: hayat
"""

# Read transposase protein list
with open(r"C:\Users\hayat\Downloads\R_files\data\mobc_plasmids.txt") as f:
    transposase_proteins = [line.strip() for line in f if line.strip()]

# Strip off the final "_<number>" to get plasmid names with transposases
transposase_plasmids = {pid.rsplit("_", 1)[0] for pid in transposase_proteins}

# Read all plasmid names
with open(r"C:\Users\hayat\Downloads\R_files\data\top_abundant_and_widespread_plasmid_names.txt") as f:
    all_plasmids = [line.strip() for line in f if line.strip()]

# Create iTOL colorstrip dataset content
output_lines = [
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tMob",
    "COLOR\t#000000",
    "LEGEND_TITLE\tMob",
    "LEGEND_SHAPES\t1\t1",
    "LEGEND_COLORS\t#351c75\t#f1c232",
    "LEGEND_LABELS\tMob+\tMob-",
    "DATA"
]

# Assign colors and labels
for plasmid in sorted(all_plasmids):
    if plasmid in transposase_plasmids:
        color = "#351c75"
        label = "Mob+"
    else:
        color = "#f1c232"
        label = "Mob-"
    output_lines.append(f"{plasmid}\t{color}\t{label}")

# Save to file
with open(r"C:\Users\hayat\Downloads\R_files\data\itol_mob.txt", "w") as out:
    out.write("\n".join(output_lines))

print("iTOL file saved as 'itol_transposase.txt'")
