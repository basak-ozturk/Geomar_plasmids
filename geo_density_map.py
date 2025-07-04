# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 10:16:59 2025

@author: hayat
"""

import pandas as pd
import plotly.express as px

# Load merged file with all plasmid presences and metadata
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmid_metadata.tsv", sep="\t")

# Drop rows with missing coordinates
df = df.dropna(subset=['Latitude', 'Longitude'])

# Group by location (round coordinates slightly to avoid jitter)
df['Latitude_round'] = df['Latitude'].round(4)
df['Longitude_round'] = df['Longitude'].round(4)

# Count unique plasmids per location
summary = df.groupby(['Latitude_round', 'Longitude_round']).agg({
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
fig.write_html("C:/Users/hayat/Downloads/R_files/graphs/plasmid_density_map.html")

print("Map saved: plasmid_density_map.html")
