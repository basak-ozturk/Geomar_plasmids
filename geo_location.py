# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:44:15 2025

@author: hayat
"""

import pandas as pd
import plotly.express as px

from geopy.distance import great_circle
import numpy as np
from sklearn.cluster import AgglomerativeClustering

def parse_lat_lon(latlon_str):
    if pd.isna(latlon_str):
        return pd.Series({'Latitude': np.nan, 'Longitude': np.nan})
    parts = latlon_str.strip().split()
    lat = float(parts[0]) * (-1 if parts[1] == 'S' else 1)
    lon = float(parts[2]) * (-1 if parts[3] == 'W' else 1)
    return pd.Series({'Latitude': lat, 'Longitude': lon})


# Load metadata
metadata = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/GEOLOCATION_HostGenus_Samples.txt", sep="\t")

# Standardize the 'Run' column (just in case)
metadata['Run'] = metadata['Run'].astype(str)

# Load plasmid list
with open("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_names.txt") as f:
    plasmids = [line.strip() for line in f if line.strip()]

# Extract Run from each plasmid
plasmid_data = []
for plasmid in plasmids:
    run = plasmid.split("_")[0]
    plasmid_data.append({"Plasmid": plasmid, "Run": run})

# Create DataFrame and merge
plasmid_df = pd.DataFrame(plasmid_data)
merged = plasmid_df.merge(metadata, on="Run", how="left")

merged = merged.rename(columns={"biome_genus": "host_genus"})

final_df = merged[["Plasmid", "Run", "lat_lon", "host_genus"]]
final_df = merged[["Plasmid", "Run", "lat_lon", "host_genus"]].copy()
final_df[['Latitude', 'Longitude']] = final_df['lat_lon'].apply(parse_lat_lon)


fig = px.scatter_geo(final_df, lat='Latitude', lon='Longitude',
                     color='host_genus', hover_name='Plasmid',
                     projection='natural earth')
fig.update_traces(marker=dict(size=10, line=dict(width=0.5, color='black')))
fig.update_layout(title="Global Occurrence Map for the Top Widespread and Abundant Plasmids", legend_title="Host Genus")
#fig.write_html("C:/Users/hayat/Downloads/R_files/graphs/plasmid_map.html")
fig.write_image("C:/Users/hayat/Downloads/R_files/graphs/plasmid_map.png", format="png", scale=2)
fig.show()



# coords = final_df[['Latitude', 'Longitude']].values

# # Calculate pairwise distance matrix (in km)
# dist_matrix = np.zeros((len(coords), len(coords)))
# for i in range(len(coords)):
#     for j in range(len(coords)):
#         dist_matrix[i,j] = great_circle(coords[i], coords[j]).km

# # Cluster plasmids by location (e.g., 3 clusters)
# cluster = AgglomerativeClustering(n_clusters=10, metric='precomputed', linkage='complete')
# labels = cluster.fit_predict(dist_matrix)

# final_df['geo_cluster'] = labels
# final_df.to_csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmid_metadata.tsv", sep="\t", index=False)
