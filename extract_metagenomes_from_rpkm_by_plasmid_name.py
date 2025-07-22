# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 15:41:00 2025

@author: hayat
"""

import pandas as pd

# Input files
plasmid_list_file = "C:/Users/hayat/Downloads/R_files/data/kristina_related_plasmids.txt"     
rpkm_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv", sep="\t", index_col=0)
output_file = "C:/Users/hayat/Downloads/R_files/data/kristina_related_plasmids_metagenomes.csv"
location_file = "C:/Users/hayat/Downloads/R_files/data/GEOLOCATION_HostGenus_Samples.txt"
# Load plasmid list
with open(plasmid_list_file, "r") as f:
    plasmid_list = [line.strip() for line in f if line.strip()]


loc_df = pd.read_csv(location_file, sep="\t")  # expects 'Run' and 'lat_lon' columns
loc_dict = loc_df.set_index("Run")["lat_lon"].to_dict()
loc_df["Run"] = loc_df["Run"].astype(str).str.strip()
rpkm_df.columns = rpkm_df.columns.str.strip()

# Group by Run and collapse all lat_lon values
loc_df_grouped = (
    loc_df.groupby("Run")["lat_lon"]
    .apply(lambda x: ";".join(sorted(set(x))))
    .reset_index()
)

# Now rebuild the lookup dictionary
loc_dict = dict(zip(loc_df_grouped["Run"], loc_df_grouped["lat_lon"]))

# Filter for valid plasmids
plasmid_list_found = [p for p in plasmid_list if p in rpkm_df.index]
missing = set(plasmid_list) - set(plasmid_list_found)
if missing:
    print(f"Warning: {len(missing)} plasmids not found in RPKM matrix: {missing}")

# Build result: plasmid, metagenomes with lat_lon
result = []
for plasmid in plasmid_list_found:
    row = rpkm_df.loc[plasmid]
    metagenomes = row[row >= 1].index.tolist()
    metagenomes_with_loc = [
        f"{mg} ({loc_dict.get(mg, 'no_location')})" for mg in metagenomes
    ]
    result.append({
        "Plasmid": plasmid,
        "Metagenomes_with_location": ";".join(metagenomes_with_loc)
    })
for mg in metagenomes:
    if mg not in loc_dict:
        print(f"Missing: {mg}")


# Save to CSV
output_df = pd.DataFrame(result)
output_df.to_csv(output_file, index=False)

print(f"Saved result to {output_file}")