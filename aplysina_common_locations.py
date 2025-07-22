# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 16:07:10 2025

@author: hayat

"""
import pandas as pd

def round_location(loc_str, decimals=2):
    # loc_str example: '45.5099 N 13.5600 E'
    try:
        parts = loc_str.split()
        lat = round(float(parts[0]), decimals)
        lat_dir = parts[1]
        lon = round(float(parts[2]), decimals)
        lon_dir = parts[3]
        return f"{lat} {lat_dir} {lon} {lon_dir}"
    except Exception as e:
        # fallback in case of unexpected format
        return loc_str

def extract_unique_locations(entry):
    locs_raw = [x.split("(", 1)[-1].rstrip(")") for x in entry.split(";")]
    locs_rounded = [round_location(loc) for loc in locs_raw]
    return sorted(set(locs_rounded))


# Load the results (assuming from CSV)
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/kristina_related_plasmids_metagenomes.csv")


df["Unique_Locations"] = df["Metagenomes_with_location"].apply(extract_unique_locations)
df["Num_Unique_Locations"] = df["Unique_Locations"].apply(len)

# Optional: print or save a simplified summary
summary_df = df[["Plasmid", "Num_Unique_Locations", "Unique_Locations"]]
summary_df.to_csv("C:/Users/hayat/Downloads/R_files/data/kristina_related_plasmid_location_summary.csv", index=False, encoding="utf-8")


