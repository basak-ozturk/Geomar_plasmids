# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 11:30:17 2025

@author: hayat
"""

import pandas as pd

# load Aplysina plasmid names
plasmid_ids = []
with open("C:/Users/hayat/Downloads/R_files/data/plasmids_in_Aplysina.tsv") as f:
    for line in f:
        first_field = line.strip().split("SRR")[0].split("ERR")[0].split("DRR")[0]
        # Try safer split (assuming first part is always the plasmid ID)
        if "_" in line:
            pid = line.strip().split("_")[0] + "_" + line.strip().split("_")[1].split("ctg")[0] + "ctg"
        else:
            continue
        plasmid_ids.append(pid)

plasmid_ids = list(set(plasmid_ids))  # Remove duplicates

# Load Kristina plasmids
with open("C:/Users/hayat/Downloads/R_files/data/ql.txt") as kf:
    kristina_plasmids = {line.strip().split("/")[-1].replace(".fasta", "") for line in kf}

# Load Run → Location table
metadata = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/GEOLOCATION_HostGenus_Samples.txt", sep="\t")  # Replace with correct filename
metadata["Run"] = metadata["Run"].astype(str).str.strip().str.upper()

# Create location labels from lat_lon
# Convert lat_lon string like "18.822567 S 147.63755 E" into numeric (rounded) coordinates

def extract_location_group(coord):
    try:
        parts = coord.strip().split()
        if len(parts) == 4:
            lat = float(parts[0]) * (-1 if parts[1] == "S" else 1)
            lon = float(parts[2]) * (-1 if parts[3] == "W" else 1)
            lat_round = round(lat, 1)
            lon_round = round(lon, 1)
            return f"{lat_round}_{lon_round}"
        else:
            return "Unknown"
    except Exception:
        return "Unknown"

# Apply function
metadata["Rounded_Location"] = metadata["lat_lon"].apply(extract_location_group)

# Create map for known locations
known_coords = sorted(set(c for c in metadata["Rounded_Location"].unique() if c != "Unknown"))
location_map = {coord: f"Location_{i+1}" for i, coord in enumerate(known_coords)}

# Assign mapped location labels
metadata["Location"] = metadata["Rounded_Location"].map(location_map)
metadata.loc[metadata["Rounded_Location"] == "Unknown", "Location"] = "Unknown_Location"



# Create a lookup dictionary
metadata_dict = metadata.set_index("Run").to_dict(orient="index")

# Build node table
rows = []
for plasmid in plasmid_ids:
    run_id = plasmid.split("_")[0].upper()
    is_kristina = plasmid in kristina_plasmids

    if is_kristina:
        location = "Kristina_Location"
        host = "Unknown"
        genus = "Unknown"
    else:
        meta = metadata_dict.get(run_id)
        if meta:
            location = meta.get("Location", "Unknown_Location")
            host = meta.get("Host", "NA")
            genus = meta.get("biome_genus", "NA")
        else:
            location = "Unknown_Location"
            host = "Unknown"
            genus = "Unknown"

    rows.append({
        "ID": plasmid,
        "Is_Kristina": is_kristina,
        "Location": location,
        "Host": host,
        "Biome_Genus": genus
    })

# Write to CSV
df_nodes = pd.DataFrame(rows)
df_nodes.to_csv("C:/Users/hayat/Downloads/R_files/data/aplysina_cytoscape_node_annotations.csv", index=False)


