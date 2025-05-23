import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors

# Complete list of allowed hosts
allowed_hosts = [
    "Agelas", "Ircinia", "Coscinoderma", "Aplysina", "Rhopaloeides", "Not_Defined", "Xestospongia",
    "Geodia", "Phyllospongia", "Spongilla", "Cliona", "Amphimedon", "Theonella", "Cinachyrella",
    "Niphates", "Haliclona", "Halichondria", "Manihinea", "Mycale", "Pseudoceratina", "Aiolochroia",
    "Clathria", "Spheciospongia", "Smenospongia", "Verongula", "Eunapius", "Acarnus", "Axinella",
    "Thoosa", "Cymbastela", "Petrosia", "Sarcotragus", "Dysidea", "Rhabdastrella"
]

# Load mapping file (Run -> Genus)
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt", sep="\t")

# Filter to only allowed hosts
df_filtered = df[df["biome_genus"].isin(allowed_hosts)]

# Create mapping Run -> Host
run_to_host = dict(zip(df_filtered["Run"], df_filtered["biome_genus"]))

# Generate distinct colors using matplotlib colormaps
cmap = plt.get_cmap("tab20b")  # use a qualitative colormap
if len(allowed_hosts) > 20:
    cmap = plt.get_cmap("tab20c")

# Custom high-contrast colors
colors = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628",
    "#f781bf", "#999999", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
    "#ffd92f", "#e5c494", "#b3b3b3", "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
    "#66a61e", "#e6ab02", "#a6761d", "#666666", "#8c564b", "#c49c94", "#fdbf6f",
    "#ccebc5", "#ffed6f", "#bc80bd", "#fb8072", "#80b1d3", "#bebada", "#ffffb3"
][:len(allowed_hosts)]


host_to_color = dict(zip(sorted(allowed_hosts), colors))

# iTOL header
itol_lines = [
    "DATASET_STYLE",
    "SEPARATOR TAB",
    "DATASET_LABEL\tHost",
    "COLOR\t#000000",
    "LEGEND_TITLE\tHost",
    "LEGEND_SHAPES\t" + "\t".join(["1"] * len(allowed_hosts)),
    "LEGEND_COLORS\t" + "\t".join([host_to_color[host] for host in sorted(allowed_hosts)]),
    "LEGEND_LABELS\t" + "\t".join(sorted(allowed_hosts)),
    "DATA"
]

# Read plasmid labels and build style lines
with open("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid_names.csv") as f:
    for line in f:
        tip = line.strip()
        prefix = tip.split("_")[0]
        host = run_to_host.get(prefix)
        if host:
            color = host_to_color[host]
            # Format: <label>    label    bold    <color>    <size>    <fontstyle>
            itol_lines.append(f"{tip}\tlabel\tbold\t{color}\t1\tnormal")

# Save
out_path = "C:/Users/hayat/Downloads/R_files/data/itol_host_annotation_STYLE.txt"
with open(out_path, "w") as out:
    out.write("\n".join(itol_lines))

print("Done! Upload 'itol_host_annotation_STYLE.txt' as a 'Style' dataset in iTOL.")
