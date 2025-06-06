# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 15:02:05 2025

@author: hayat
"""

import pandas as pd

# Load main dataframe
df = pd.read_csv(r"C:\Users\hayat\Downloads\R_files\data\widespread_with_integrase_and_length.csv")

# Load Run-to-host-genus mapping
accession_df = pd.read_csv(r"C:\Users\hayat\Downloads\R_files\data\Accessions_to_SpongeGenusIDs.txt", sep="\t")

# Build dictionary: Run → biome_genus
accession_to_host = dict(zip(accession_df["Run"], accession_df["biome_genus"]))

# Extract metagenome accession from plasmid name (before first underscore)
df["Run"] = df["Plasmid"].str.extract(r"(.+?)_")[0]

# Map to host genus
df["Host"] = df["Run"].map(accession_to_host)

# Define host-to-color map
host_colors = {
    "Acarnus": "#1f77b4", "Agelas": "#ff7f0e", "Aiolochroia": "#2ca02c", "Amphimedon": "#d62728",
    "Aplysina": "#9467bd", "Axinella": "#8c564b", "Baikalospongia": "#e377c2", "Chondrilla": "#7f7f7f",
    "Cinachyrella": "#bcbd22", "Clathria": "#17becf", "Cliona": "#aec7e8", "Coscinoderma": "#ffbb78",
    "Crella": "#98df8a", "Cymbastela": "#ff9896", "Dysidea": "#c5b0d5", "Ephydatia": "#c49c94",
    "Eunapius": "#f7b6d2", "Geodia": "#dbdb8d", "Halichondria": "#9edae5", "Haliclona": "#393b79",
    "Haplosclerida": "#637939", "Hymedesmia": "#8c6d31", "Ircinia": "#FFFF00",  
    "Isodictya": "#843c39", "Lamellodysidea": "#7b4173", "Lanthella": "#5254a3", "Leucetta": "#6b6ecf",
    "Lophophysema": "#9c9ede", "Manihinea": "#637939", "Mycale": "#8ca252", "Myxilla": "#b5cf6b",
    "Niphates": "#cedb9c", "Not_Defined": "#bd9e39", "Pericharax": "#e7ba52", "Petrosia": "#e7969c",
    "Phyllospongia": "#d6616b", "Pseudoceratina": "#ad494a", "Rhabdastrella": "#843c39",
    "Rhopaloeides": "#d94801", "Sarcotragus": "#f16913", "Scopalina": "#fd8d3c", "Smenospongia": "#fdae6b",
    "Spheciospongia": "#fdd0a2", "Spongilla": "#31a354", "Stylissa": "#74c476", "Tedaniidae": "#a1d99b",
    "Theonella": "#c7e9c0", "Thoosa": "#6baed6", "Verongula": "#9ecae1", "Xestospongia": "#c6dbef"
}


# Filter host_colors to include only those hosts in the data
hosts_in_data = set(df["Host"].unique())
filtered_host_colors = {host: color for host, color in host_colors.items() if host in hosts_in_data}

# Assign colors based on filtered dictionary
df["Host_Color"] = df["Host"].map(filtered_host_colors).fillna("#aaaaaa")  # fallback grey
df["Host"] = df["Host"].fillna("Unknown")

# Save combined data
df.to_csv("C:/Users/hayat/Downloads/R_files/data/combined_plasmid_host_integrase.csv", index=False)

# Write iTOL symbol file
with open("C:/Users/hayat/Downloads/R_files/data/itol_host_symbols.txt", "w") as f:
    f.write("DATASET_SYMBOL\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tHost\n")
    f.write("COLOR\t#000000\n")
    f.write("LEGEND_TITLE\tHost\n")
    f.write("LEGEND_SHAPES\t" + "\t".join(["2"] * len(filtered_host_colors)) + "\n")
    f.write("LEGEND_COLORS\t" + "\t".join(filtered_host_colors.values()) + "\n")
    f.write("LEGEND_LABELS\t" + "\t".join(filtered_host_colors.keys()) + "\n")
    f.write("DATA\n")
    for _, row in df.iterrows():
        plasmid_id = row["Plasmid"]
        color = row["Host_Color"]
        label = row["Host"]
        f.write(f"{plasmid_id}\t2\t10\t{color}\t1\t-1\t{label}\n")



# Build integrase colorstrip
# Build plasmid length bar plot
with open("C:/Users/hayat/Downloads/R_files/data/itol_plasmid_length.txt", "w") as f:
    f.write("DATASET_SIMPLEBAR\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tPlasmid Length\n")
    f.write("COLOR\t#808080\n")
    f.write("WIDTH\t25\n")
    f.write("BORDER_WIDTH\t1\n")
    f.write("BORDER_COLOR\t#000000\n")
    f.write("DATA\n")
    for _, row in df.iterrows():
        f.write(f"{row['Plasmid']}\t{row['Sequence_Length']}\n")


