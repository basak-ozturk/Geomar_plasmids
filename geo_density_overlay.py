import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
#import plotly.io as pio


# === 1. Load density summary ===
summary = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_density_summary.csv")

# === 2. Load target plasmid locations ===
loc_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/SRR10694904_497ctg_locations.tsv", sep="\t")
loc_df = loc_df.dropna(subset=["Latitude", "Longitude"])

# === 3. Base density layer ===
density_trace = go.Scattergeo(
    lat=summary["Latitude_round"],
    lon=summary["Longitude_round"],
    mode="markers",
    marker=dict(
        size=summary["Plasmid_Count"],
        color=summary["Plasmid_Count"],
        colorscale='Greys',
        cmin=summary["Plasmid_Count"].min(),
        cmax=summary["Plasmid_Count"].max(),
        colorbar=dict(title="Plasmid Count"),
        line=dict(width=1, color='red'),
        sizemode='area',
        sizeref=2.*summary["Plasmid_Count"].max()/25**2,  # mimic px.scatter_geo size_max=30
        sizemin=2
    ),
    name="Plasmid Density",
    hovertemplate="Plasmid Count: %{marker.size}<extra></extra>"
)

# === 4. Overlay SRR10694904_497ctg points ===
host_colors = px.colors.qualitative.Dark24
unique_hosts = loc_df["host_genus"].dropna().unique()
color_map = {host: host_colors[i % len(host_colors)] for i, host in enumerate(unique_hosts)}
point_colors = loc_df["host_genus"].map(color_map)

highlight_traces = []

for i, host in enumerate(unique_hosts):
    host_df = loc_df[loc_df["host_genus"] == host]
    trace = go.Scattergeo(
        lat=host_df["Latitude"],
        lon=host_df["Longitude"],
        mode="markers",
        marker=dict(
            size=5,
            color=color_map[host],
            line=dict(width=1, color="black"),
            symbol="circle"
        ),
        name=host,  # <- this ensures legend entry per host
        hovertemplate="Host Genus: %{text}<extra></extra>",
        text=host_df["host_genus"]
    )
    highlight_traces.append(trace)


# === 5. Combine and style ===
fig = go.Figure(data=[density_trace, *highlight_traces])

fig.update_layout(
    title="Plasmid Density + SRR10694904_497ctg Overlay",
    geo=dict(
        showland=True,
        landcolor="rgb(243, 243, 243)",
        projection_type="natural earth"
    ),
    margin=dict(r=0, t=40, l=0, b=0)
)

# === 6. Export ===
fig.write_image("C:/Users/hayat/Downloads/R_files/graphs/SRR10694904_overlay_map.svg")

fig.write_html("C:/Users/hayat/Downloads/R_files/graphs/SRR10694904_overlay_matched_style.html")
print("Overlay map saved with matching style.")
