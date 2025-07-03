# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 12:15:20 2025

@author: hayat
"""

def get_plasmid_name(protein_id):
    # Remove the last underscore and trailing number
    # Example: ERR5004223_7326ctg_1 -> ERR5004223_7326ctg
    return "_".join(protein_id.split("_")[:-1])

# Load mobc proteins and get unique plasmids
with open("C:/Users/hayat/Downloads/R_files/data/widespread_plasmids_MobC_proteins.txt") as f:
    mobc_proteins = [line.strip() for line in f if line.strip()]

mobc_plasmids = set(get_plasmid_name(p) for p in mobc_proteins)

# Load oriT plasmids
with open("C:/Users/hayat/Downloads/R_files/data/oriT_alloriT_blast_results_names.txt") as f:
    orit_plasmids = set(line.strip() for line in f if line.strip())

# Load all plasmids
with open("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_names.txt") as f:
    all_plasmids = set(line.strip() for line in f if line.strip())

# Compare mobc plasmids to oriT plasmids (intersection)
mobc_and_orit = mobc_plasmids.intersection(orit_plasmids)

# Compare all plasmids to oriT plasmids
all_and_orit = all_plasmids.intersection(orit_plasmids)

# Print results
print("Unique plasmids with mobc:", len(mobc_plasmids))
print("Unique plasmids with oriT:", len(orit_plasmids))
print("Mobc plasmids that also have oriT:", len(mobc_and_orit))
print("All plasmids that also have oriT:", len(all_and_orit))

print("\nMobc plasmids with oriT:")
for p in sorted(mobc_and_orit):
    print(p)

print("\nAll plasmids with oriT:")
for p in sorted(all_and_orit):
    print(p)

# Save mobc plasmids to file
with open("mobc_plasmids.txt", "w") as out_f:
    for plasmid in sorted(mobc_plasmids):
        out_f.write(plasmid + "\n")

print(f"Saved {len(mobc_plasmids)} unique mobc plasmids to mobc_plasmids.txt")