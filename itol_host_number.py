import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load main dataframe (tip IDs in first column)
df = pd.read_csv(r"C:\Users\hayat\Downloads\R_files\data\widespread_with_integrase_and_length.csv")

# Extract tip IDs from the first column (assumed to be named or use .iloc)
tree_ids = set(df.iloc[:, 0].astype(str).tolist())

# Load plasmid-host heatmap
heatmap_df = pd.read_csv(r"C:\Users\hayat\Downloads\R_files\data\plasmids_vs_sponge_genera.csv")

# Prepare host genera list (exclude 'Genome')
host_genera = list(heatmap_df.columns)
host_genera.remove("Genome")

max_hosts = len(host_genera)
# Prepare color map from light gray to blue
cmap = plt.get_cmap("Blues")
min_hosts = 0
max_hosts = len(host_genera)

def get_color(host_count):
    # Normalize host count between 0 and 1
    norm = (host_count - min_hosts) / (max_hosts - min_hosts) if max_hosts > 0 else 0
    rgba = cmap(norm)
    # Convert RGBA to hex
    return mcolors.to_hex(rgba)

min_hosts_found = heatmap_df[host_genera].gt(0).sum(axis=1).min()
max_hosts_found = heatmap_df[host_genera].gt(0).sum(axis=1).max()


with open("C:/Users/hayat/Downloads/R_files/data/itol_host_count_colorstrip.txt", "w") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tHost_Count\n")
    f.write("COLOR\t#000000\n")  # Border color of strip
    f.write("LEGEND_TITLE\tNumber of Hosts\n")
    f.write("LEGEND_SHAPES\t1\t1\n")
    f.write("LEGEND_COLORS\t#d3d3d3\t#1f77b4\n")
    f.write(f"LEGEND_LABELS\t{min_hosts_found}\t{max_hosts_found}\n")
    f.write("DATA\n")

    for _, row in heatmap_df.iterrows():
        plasmid_id = row["Genome"]
        if plasmid_id not in tree_ids:
            # Skip plasmid IDs not in the tree
            continue
        host_count = sum(row[host] > 0 for host in host_genera)
        if host_count == 0:
            continue
        color = get_color(host_count)
        f.write(f"{plasmid_id}\t{color}\n")

# Save list of plasmids with number of hosts
plasmid_host_counts = []

for _, row in heatmap_df.iterrows():
    plasmid_id = row["Genome"]
    if plasmid_id not in tree_ids:
        continue
    host_count = sum(row[host] > 0 for host in host_genera)
    if host_count == 0:
        continue
    plasmid_host_counts.append((plasmid_id, host_count))

# Create and save DataFrame
host_count_df = pd.DataFrame(plasmid_host_counts, columns=["Plasmid", "Host_Count"])
host_count_df.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_host_counts.csv", index=False)

# Load host count file
host_count_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_host_counts.csv")

# Load integrase and length file
info_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/widespread_with_integrase_and_length.csv")

# Merge on 'Plasmid'
merged_df = pd.merge(info_df, host_count_df, on="Plasmid", how="left")

# Fill missing Host_Count with 0 (for plasmids not found in host count file)
merged_df["Host_Count"] = merged_df["Host_Count"].fillna(0).astype(int)

# Save to new file
merged_df.to_csv("C:/Users/hayat/Downloads/R_files/data/combined_plasmid_info_with_host_counts.csv", index=False)
