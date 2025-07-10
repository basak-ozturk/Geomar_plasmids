import pandas as pd


input_tsv = "C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmid_metadata.tsv"
output_itol = "C:/Users/hayat/Downloads/R_files/data/itol_geo_cluster_symbols_extended.txt"

# Load metadata
df = pd.read_csv(input_tsv, sep="\t")

if 'geo_cluster' not in df.columns:
    raise ValueError("Missing 'geo_cluster' column")

clusters = sorted(df['geo_cluster'].unique())
n_clusters = len(clusters)

# Color Universal Design (colorblind-friendly) palette
cud_colors = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#000000", "#999933"
]

colors = [cud_colors[i % len(cud_colors)] for i in range(n_clusters)]
cluster_to_color = {c: colors[i] for i, c in enumerate(clusters)}

# Parameters
shape_code = 4  # Star
size = 10
fill = 1        # 1 = filled
position = 1    # Symbol position

# Legend entries
legend_labels = [f"Cluster {c}" for c in clusters]
legend_shapes = [str(shape_code)] * n_clusters
legend_colors = colors

# Write iTOL file
with open(output_itol, "w", encoding="utf-8", newline='') as f:
    f.write("DATASET_SYMBOL\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tGeo Cluster (Extended)\n")
    f.write("COLOR\t#000000\n")
    f.write("LEGEND_TITLE\tGeo Cluster\n")
    f.write(f"LEGEND_SHAPES\t{'\t'.join(legend_shapes)}\n")
    f.write(f"LEGEND_COLORS\t{'\t'.join(legend_colors)}\n")
    f.write(f"LEGEND_LABELS\t{'\t'.join(legend_labels)}\n")
    f.write("DATA\n")

    for _, row in df.iterrows():
        plasmid = row['Plasmid']
        cluster = row['geo_cluster']
        color = cluster_to_color[cluster]
        label = f"Cluster {cluster}"
        # Write: ID, symbol, size, color, fill (1), position, label
        f.write(f"{plasmid}\t{shape_code}\t{size}\t{color}\t{fill}\t{position}\t{label}\n")

print(f"✅ iTOL symbol dataset saved: {output_itol}")
