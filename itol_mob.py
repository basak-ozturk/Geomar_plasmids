# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 09:02:08 2025

@author: hayat
"""

with open("C:/Users/hayat/Downloads/R_files/data/widespread_plasmids_relaxase_proteins.txt") as f:
    mob_ids = [line.strip().rsplit('_', 1)[0] for line in f]

mob_plasmids = set(mob_ids)


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
    if plasmid in mob_plasmids:
        color = "#351c75"
        label = "Mob+"
    else:
        color = "#f1c232"
        label = "Mob-"
    output_lines.append(f"{plasmid}\t{color}\t{label}")

# Save to file
with open(r"C:\Users\hayat\Downloads\R_files\data\itol_mob.txt", "w") as out:
    out.write("\n".join(output_lines))

print("iTOL file saved as 'itol_mob.txt'")
