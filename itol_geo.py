import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

input_tsv = "C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmid_metadata.tsv"
output_itol = "C:/Users/hayat/Downloads/R_files/data/itol_geo_cluster_symbols_extended.txt"

df = pd.read_csv(input_tsv, sep="\t")

if 'geo_cluster' not in df.columns:
    raise ValueError("Missing 'geo_cluster' column")

clusters = sorted(df['geo_cluster'].unique())
n_clusters = len(clusters)

cmap = plt.get_cmap('tab10')
colors = [mcolors.rgb2hex(cmap(i % 10)) for i in range(n_clusters)]
cluster_to_color = {c: colors[i] for i, c in enumerate(clusters)}

shape_code = 4  
size = 10
border_color = "#000000"
border_width = 1

legend_labels = [f"Cluster {c}" for c in clusters]
legend_shapes = [str(shape_code)] * n_clusters
legend_colors = colors

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
        f.write(f"{plasmid}\t{shape_code}\t{size}\t{color}\t{border_color}\t{border_width}\t{label}\n")

print(f"✅ Extended iTOL DATASET_SYMBOL file written to: {output_itol}")
