import pandas as pd
import plotly.express as px
import glob
import os

# Define folders
data_folder = "C:/Users/hayat/Downloads/R_files/data/"
output_folder = "C:/Users/hayat/Downloads/R_files/graphs/"
os.makedirs(output_folder, exist_ok=True)

# Get all *_locations.tsv files
tsv_files = glob.glob(os.path.join(data_folder, "*_locations.tsv"))

# Loop through each file
for file_path in tsv_files:
    # Extract plasmid name from filename
    filename = os.path.basename(file_path)
    plasmid_name = filename.replace("_locations.tsv", "")
    
    # Load the metadata
    metadata = pd.read_csv(file_path, sep="\t")
    
    # Drop rows with missing coordinates
    metadata = metadata.dropna(subset=['Latitude', 'Longitude'])

    # If there's nothing left to plot, skip
    if metadata.empty:
        print(f"Skipping {filename} — all rows missing coordinates.")
        continue

    # Plot
    fig = px.scatter_geo(metadata,
                         lat='Latitude',
                         lon='Longitude',
                         color='host_genus',
                         hover_name='Plasmid',
                         projection='natural earth')

    # Save HTML
    fig.update_traces(marker=dict(size=10, line=dict(width=0.5, color='black')))
    fig.update_layout(title=f"{filename} — Occurrence Map",
                      legend_title="Host Genus")
    output_file = os.path.join(output_folder, f"{plasmid_name}_plasmid_map.html")
    fig.write_html(output_file)
    print(f"Saved map for {plasmid_name} to {output_file}")
    # Save PNG
    output_file_png = os.path.join(output_folder, f"{plasmid_name}_plasmid_map.png")
    fig.write_image(output_file_png, format="png", scale=2)  # scale=2 makes it higher resolution
    print(f"Saved PNG map for {plasmid_name} to {output_file_png}")


# Combine all files into one DataFrame
all_data = []
output_file2 = "C:/Users/hayat/Downloads/R_files/graphs/combined_plasmid_map.html"
for file_path in tsv_files:
    df = pd.read_csv(file_path, sep="\t")
    plasmid_name = os.path.basename(file_path).replace("_locations.tsv", "")
    df['Plasmid'] = plasmid_name  # Ensure consistent plasmid label
    df = df.dropna(subset=['Latitude', 'Longitude'])
    all_data.append(df)

# Concatenate into single DataFrame
combined_df = pd.concat(all_data, ignore_index=True)

# Plot with different shapes per plasmid and larger markers
fig = px.scatter_geo(combined_df,
                     lat='Latitude',
                     lon='Longitude',
                     color='host_genus',
                     symbol='Plasmid',  # Unique symbol per plasmid
                     hover_name='Plasmid',
                     projection='natural earth',
                     size_max=20)

# Increase symbol size manually (Plotly doesn't scale size well for geo)
fig.update_traces(marker=dict(size=10, line=dict(width=0.5, color='black')))

# Adjust layout
fig.update_layout(title="Combined Plasmid Occurrence Map",
                  legend_title="Host Genus / Plasmid")

# Save to HTML
fig.write_html(output_file2)
#print(f"Combined map saved to {output_file2}")


output_file2_png = "C:/Users/hayat/Downloads/R_files/graphs/combined_plasmid_map.png"
fig.write_image(output_file2_png, format="png", scale=2)
print(f"Combined PNG map saved to {output_file2_png}")
