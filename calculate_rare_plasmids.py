import pandas as pd

# Load data
df = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv",
    sep="\t",
    index_col=0
)

# Load plasmid list
with open("C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids_names.txt", "r") as f:
    plasmid_names_list = [line.strip() for line in f if line.strip()]

# Keep only plasmids of interest
df = df.loc[df.index.isin(plasmid_names_list)]

# Presence/absence (RPKM ≥ 1 counts as present)
presence = df >= 1

# Count in how many metagenomes each plasmid is present
plasmid_presence_counts = presence.sum(axis=1)

# Calculate fraction of metagenomes for each plasmid
fraction_present = plasmid_presence_counts / presence.shape[1]

# Find plasmids present in ≤1% of metagenomes
rare_plasmids = fraction_present[fraction_present <= 0.01]

print(f"Number of plasmids present in ≤1% of metagenomes: {rare_plasmids.shape[0]}")
