# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 10:58:48 2025

@author: hayat
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:44:15 2025

@author: hayat
"""

import pandas as pd
import plotly.express as px
import plotly.colors as pc
#from geopy.distance import great_circle
import numpy as np
#from sklearn.cluster import AgglomerativeClustering

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
with open("C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids_names.txt") as f:
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
final_df[['Latitude', 'Longitude']] = final_df['lat_lon'].apply(parse_lat_lon)

# Drop rows with missing coordinates before plotting or clustering
final_df = final_df.dropna(subset=['Latitude', 'Longitude'])

color_sequence = pc.qualitative.Alphabet

fig = px.scatter_geo(final_df, lat='Latitude', lon='Longitude',
                     color='host_genus', hover_name='Plasmid',
                     projection='natural earth',
                     color_discrete_sequence=color_sequence)
fig.update_traces(marker=dict(size=10, line=dict(width=0.5, color='black')))

fig.write_html("C:/Users/hayat/Downloads/R_files/graphs/all_plasmid_map.html")


final_df['Latitude_round'] = final_df['Latitude'].round(4)
final_df['Longitude_round'] = final_df['Longitude'].round(4)

# Count unique plasmids per location
summary = final_df.groupby(['Latitude_round', 'Longitude_round']).agg({
    'Plasmid': pd.Series.nunique
}).reset_index().rename(columns={'Plasmid': 'Plasmid_Count'})

# Plot
fig = px.scatter_geo(summary,
                     lat='Latitude_round',
                     lon='Longitude_round',
                     size='Plasmid_Count',
                     color='Plasmid_Count',
                     color_continuous_scale='Viridis',
                     projection='natural earth',
                     title='Unique Plasmids per Sampling Location')

fig.update_traces(marker=dict(line=dict(width=0)))
fig.update_layout(geo=dict(showland=True, landcolor="rgb(243, 243, 243)"),
                  margin={"r":0,"t":40,"l":0,"b":0})
fig.update_layout(legend_title="Host Genus / Plasmid")
# Save HTML
fig.write_html("C:/Users/hayat/Downloads/R_files/graphs/all_plasmid_density_map.html")

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
final_df.to_csv("C:/Users/hayat/Downloads/R_files/data/all_plasmid_metadata.tsv", sep="\t", index=False)
