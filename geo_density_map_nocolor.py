import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# Load merged file with all plasmid presences and metadata
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmid_metadata.tsv", sep="\t")
df = df.dropna(subset=['Latitude', 'Longitude'])

# Round coordinates slightly to avoid jitter
df['Latitude_round'] = df['Latitude'].round(4)
df['Longitude_round'] = df['Longitude'].round(4)

# Count unique plasmids per location
summary = df.groupby(['Latitude_round', 'Longitude_round']).agg({
    'Plasmid': pd.Series.nunique
}).reset_index().rename(columns={'Plasmid': 'Plasmid_Count'})

# Scale size for visibility
summary['Scaled_Count'] = summary['Plasmid_Count'] * 3  # adjust this factor as needed

# Base map using plotly express
fig = px.scatter_geo(summary,
                     lat='Latitude_round',
                     lon='Longitude_round',
                     size='Scaled_Count',
                     projection='natural earth',
                     title='Unique Plasmids per Sampling Location',
                     size_max=30)

# Convert to figure and update marker styling
fig.update_traces(
    marker=dict(
        color='red',
        line=dict(color='red', width=1),
        symbol='circle-open'
    ),
    showlegend=False  # hide auto legend
)

# --- Custom size legend ---
# Define some example plasmid counts to show in legend
legend_sizes = [5, 10, 20]
legend_scale_factor = 3

# Add one fake trace per legend item
for count in legend_sizes:
    fig.add_trace(go.Scattergeo(
        lon=[None], lat=[None],  # no visible point
        mode='markers',
        marker=dict(
            size=count * legend_scale_factor,
            color='red',
            line=dict(color='red', width=1),
            symbol='circle-open'
        ),
        showlegend=True,
        name=f'{count} plasmids'
    ))

# Layout tweaks
fig.update_layout(
    geo=dict(showland=True, landcolor="rgb(243, 243, 243)"),
    margin={"r":0, "t":40, "l":0, "b":0},
    legend_title_text='Density at Location'
)

# Save
fig.write_html("C:/Users/hayat/Downloads/R_files/graphs/plasmid_density_map_with_custom_legend.html")
print("Map saved: plasmid_density_map_with_custom_legend.html")
