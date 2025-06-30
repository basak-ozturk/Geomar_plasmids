import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.decomposition import PCA
#from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu

# Load your data
df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/CoverM_MAPPING_rpkm_Plasmid_Contigs_ouput.tsv", sep="\t", index_col=0)


# Load the plasmid names from the text file into a Python list
with open("C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids_names.txt", "r") as f:
    plasmid_names_list = [line.strip() for line in f if line.strip()]
df = df.loc[df.index.isin(plasmid_names_list)]



# Compute each plasmid’s maximum (log-transformed) RPKM across all metagenomes
max_rpkm = df.max(axis=1)

# Inspect the distribution of those maxima
plt.figure(figsize=(6,4))
max_rpkm.hist(bins=50)
plt.title("Distribution of max(log2-RPKM+1) per plasmid")
plt.xlabel("max(log2-RPKM+1)")
plt.ylabel("count of plasmids")
plt.tight_layout()
plt.show()

print(df.describe())

plt.show()

nonzero_values = df.values[df.values > 0]
log_values = np.log10(nonzero_values + 1)

plt.figure(figsize=(8, 4))
plt.hist(log_values, bins=100, color='steelblue', edgecolor='k')
plt.axvline(np.log10(1+ 1), color='red', linestyle='--', label='RPKM = 1')
plt.xlabel("log10(RPKM + 1)")
plt.ylabel("Frequency")
plt.title("Distribution of non-zero RPKM values")
plt.legend()

# Set x-axis ticks every 0.1 between 0 and max(log_values)
xticks = np.arange(0, np.ceil(log_values.max()) + 0.5, 0.5)
plt.xticks(xticks)

plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_RPKM_log10.png", dpi=300, bbox_inches="tight")

plt.show()

min_presence = int(0.01 * df.shape[1])  #optional filtering to get rid of some of the data sparsity
presence = df >= 1
filtered_df = df[presence.sum(axis=1) >= min_presence]
# Drop metagenomes that have only zeros after plasmid filtering
filtered_df = filtered_df.loc[:, (filtered_df > 0).any(axis=0)]

print(filtered_df.shape)

# RPKM distribution per metagenome

# Reset index to bring plasmid names into a column
filtered_df_reset = filtered_df.reset_index().rename(columns={'index': 'Plasmid'})

# Melt so each row is a plasmid's RPKM in a metagenome
melted_df = filtered_df_reset.melt(id_vars='Genome', var_name='Metagenome', value_name='RPKM')

# Optional: drop zero-RPKM values for better clarity in plots
melted_df_nonzero = melted_df[melted_df["RPKM"] > 0]

# Boxplot of plasmid abundance per metagenome
plt.figure(figsize=(14, 6))
sns.boxplot(data=melted_df_nonzero, x="Metagenome", y="RPKM")
plt.xticks(rotation=90)
plt.yscale("log")
plt.title("Distribution of Plasmid RPKM Values per Metagenome")
plt.tight_layout()
plt.show()

# Plot top N ranked RPKMs for each metagenome
plt.figure(figsize=(10, 6))

for sample in filtered_df.index:
    sorted_rpkm = np.sort(filtered_df.loc[sample])[::-1]  # sort descending
    plt.plot(range(1, len(sorted_rpkm)+1), sorted_rpkm, label=sample)

plt.yscale('log')
plt.xlabel('Plasmid Rank (1 = most abundant)')
plt.ylabel('RPKM')
plt.title('Ranked Plasmid Abundance per Metagenome')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

binary_matrix = (filtered_df >= 1).astype(int)

# PCA on plasmids (rows)
pca = PCA(n_components=2)
coords = pca.fit_transform(binary_matrix.values)

plt.figure(figsize=(8, 6))
plt.scatter(coords[:, 0], coords[:, 1], alpha=0.5, s=10)
plt.title("PCA of plasmid presence across metagenomes")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.show()

explained_variance = pca.explained_variance_ratio_
for i, var in enumerate(explained_variance, 1):
    print(f"PC{i} explains {var:.2%} of the variance")

kmeans = KMeans(n_clusters=2, random_state=42)
clusters = kmeans.fit_predict(coords)  # coords = pca-transformed data

plt.figure(figsize=(8, 6))
scatter = plt.scatter(coords[:, 0], coords[:, 1], c=clusters, cmap='Set1', s=20, alpha=0.7)

plt.xlabel(f"PC1 ({explained_variance[0]:.1%} variance explained)")
plt.ylabel(f"PC2 ({explained_variance[1]:.1%} variance explained)")
plt.title("KMeans Clustering on PCA of Plasmid Presence Across Metagenomes")
#plt.colorbar(scatter, label='Cluster')
plt.tight_layout()

#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/PCA_plasmid_abundance_final_plasmids.png", dpi=300, bbox_inches="tight")

plt.show()

# Add cluster labels to the binary matrix
clustered_df = binary_matrix.copy()
clustered_df["cluster"] = clusters

clustered_df["metagenome_presence_count"] = clustered_df.drop(columns="cluster").sum(axis=1)

sns.boxplot(x="cluster", y="metagenome_presence_count", data=clustered_df.reset_index())
plt.title("Plasmid Presence Across Metagenomes by Cluster")
plt.ylabel("Number of Metagenomes where Plasmid is Present")
plt.xlabel("KMeans Cluster")
plt.show()

summary = clustered_df.groupby("cluster")["metagenome_presence_count"].describe()
print(summary)

mean_profiles = clustered_df.groupby("cluster").mean(numeric_only=True).drop(columns=["metagenome_presence_count"])

plt.figure(figsize=(14, 6))  # Set the figure size here
sns.heatmap(mean_profiles.T, cmap="viridis")
plt.title("Average presence of each cluster across metagenomes")
plt.xlabel("Cluster")
plt.ylabel("Metagenome")
plt.yticks([])

plt.tight_layout()

#plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_presence_across_metagenomes_final_plasmids.png", dpi=300, bbox_inches="tight")

plt.show()

# Add a column that counts presence across metagenomes
clustered_df["metagenome_presence_count"] = (clustered_df.drop(columns=["cluster"]) > 0).sum(axis=1)

# Compare mean presence counts across clusters
print(clustered_df.groupby("cluster")["metagenome_presence_count"].describe())

top_global = clustered_df.sort_values("metagenome_presence_count", ascending=False).head(20)
print(top_global)

top_global.to_csv("C:/Users/hayat/Downloads/R_files/data/top_global_plasmids_final_plasmids.csv")

most_abundant_by_sample = filtered_df.idxmax()
abundance_values = filtered_df.max()

abundance_summary = pd.DataFrame({
    "most_abundant_plasmid": most_abundant_by_sample,
    "RPKM": abundance_values
})
print(abundance_summary.head())

# Reset index and rename the index column
abundance_summary_reset = abundance_summary.reset_index()
abundance_summary_reset = abundance_summary_reset.rename(columns={"index": "Metagenome"})

# Save with headers
abundance_summary_reset.to_csv("C:/Users/hayat/Downloads/R_files/data/most_abundant_plasmid_per_metagenome_final_plasmids.csv", index=False)


# Load the host metadata file
host_info = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/Accessions_to_SpongeGenusIDs.txt", sep="\t")


# Merge with the abundance summary
merged_df = abundance_summary_reset.merge(host_info, left_on="Metagenome", right_on="Run", how="left")

# Drop the redundant 'Run' column if desired
merged_df = merged_df.drop(columns=["Run"])
#merged_df["biome_genus"] = merged_df["biome_genus"].fillna("Spongilla")

# Save the merged result
merged_df.to_csv("C:/Users/hayat/Downloads/R_files/data/abundance_with_host_final_plasmids.csv", index=False)

print(merged_df.head())


# Merge host info into the transposed filtered_df
filtered_with_genus = filtered_df.T.merge(
    host_info, left_index=True, right_on="Run", how="left"
)

# Optional: drop "Run" column
filtered_with_genus = filtered_with_genus.drop(columns=["Run"])

# Number of metagenomes per genus

# Count how many metagenomes belong to each genus
genus_counts = host_info["biome_genus"].value_counts().reset_index()
genus_counts.columns = ["Sponge Genus", "Number of Metagenomes"]

# Sort by count for better visuals
genus_counts = genus_counts.sort_values("Number of Metagenomes", ascending=False)

# Plot
plt.figure(figsize=(10, 6))
sns.barplot(data=genus_counts, x="Sponge Genus", y="Number of Metagenomes", palette="viridis")

plt.title("Number of Metagenomes per Sponge Genus")
plt.ylabel("Number of Metagenomes")
plt.xlabel("Sponge Genus")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()


plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/metagenome_count_per_genus_final_plasmids.png", dpi=300, bbox_inches="tight")

plt.show()

# Number of plasmids per genus

# Transpose filtered_df to map metagenomes to plasmids
plasmid_presence = (filtered_df >= 1).astype(int)  # binary presence/absence
plasmid_presence_T = plasmid_presence.T  # Now rows = metagenomes

#  Merge presence info with genus info
plasmid_with_genus = plasmid_presence_T.merge(host_info, left_index=True, right_on="Run", how="left")

# Group by genus and count plasmid presence
# For each genus, get the plasmids present in any of its metagenomes
plasmids_per_genus = (
    plasmid_with_genus
    .groupby("biome_genus")
    .apply(lambda df: (df.drop(columns=["Run", "biome_genus"]).sum(axis=0) > 0).sum())
    .reset_index()
)

plasmids_per_genus.columns = ["Sponge Genus", "Number of Plasmids"]
plasmids_per_genus = plasmids_per_genus.sort_values("Number of Plasmids", ascending=False)


plt.figure(figsize=(10, 6))
sns.barplot(data=plasmids_per_genus, x="Sponge Genus", y="Number of Plasmids", palette="mako")

plt.title("Number of Plasmids per Sponge Genus")
plt.ylabel("Number of Plasmids (RPKM ≥ 1)")
plt.xlabel("Sponge Genus")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_count_per_genus_final_plasmids.png", dpi=300, bbox_inches="tight")

plt.show()

# Plasmid abundance (RPKM) per genus


# total or mean RPKM per genus
genus_total_rpkm = merged_df.groupby("biome_genus")["RPKM"].sum().sort_values(ascending=False)


log_total_rpkm = np.log10(genus_total_rpkm)

plt.figure(figsize=(12, 6))
sns.barplot(x=log_total_rpkm.index, y=genus_total_rpkm.values, palette="muted")
plt.xticks(rotation=45, ha="right")
plt.ylabel("log10(Total RPKM)")
plt.title("Total Plasmid Abundance per Sponge Genus")
plt.tight_layout()

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/rpkm_count_per_genus_final_plasmids.png", dpi=300, bbox_inches="tight")

plt.show()

plt.figure(figsize=(12, 6))
sns.barplot(x=log_total_rpkm.index, y=log_total_rpkm.values, palette="muted")
plt.xticks(rotation=45, ha="right")
plt.ylabel("log10(Total RPKM)")
plt.title("Log-Transformed Total Plasmid Abundance per Sponge Genus")
plt.tight_layout()

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/rpkm_count_per_genus_log_final_plasmids.png", dpi=300, bbox_inches="tight")

plt.show()

#  Compute total or mean RPKM per metagenome
merged_df["Total_RPKM"] = merged_df["RPKM"]  # Already computed, but we can rename for clarity

plt.figure(figsize=(10, 6))

# Scatter with jittered x-values per genus
sns.stripplot(data=merged_df, x="biome_genus", y="Total_RPKM", jitter=True, hue="biome_genus", palette="tab10", dodge=False)

plt.yscale("log")  # Optional, if you have high variance
plt.xlabel("Sponge Genus")
plt.ylabel("Total RPKM per Metagenome")
plt.title("Plasmid Abundance per Metagenome, Colored by Host Genus")
plt.legend(title="Genus", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/scatter_metagenome_rpkm_by_genus_final_plasmids.png", dpi=300, bbox_inches="tight")
plt.show()

# Create boxplot

# Filter to genera with at least 3 metagenomes
filtered_genus_counts = genus_counts[genus_counts["Number of Metagenomes"] >= 3].copy()

# Create label
filtered_genus_counts["Label"] = filtered_genus_counts["Sponge Genus"] + " (n=" + filtered_genus_counts["Number of Metagenomes"].astype(str) + ")"

# Map label to merged_df
genus_to_label = dict(zip(filtered_genus_counts["Sponge Genus"], filtered_genus_counts["Label"]))
merged_df["genus_with_count"] = merged_df["biome_genus"].map(genus_to_label)

# Filter merged_df to include only those genera
filtered_merged_df = merged_df[merged_df["genus_with_count"].notna()]

# Define x-axis order
category_order = filtered_genus_counts["Label"].tolist()

# Plot
plt.figure(figsize=(12, 6))
sns.boxplot(
    data=filtered_merged_df,
    x="genus_with_count",
    y="Total_RPKM",
    #palette="tab10",
    order=category_order
)

plt.yscale("log")
plt.xlabel("Sponge Genus")
plt.ylabel("Total RPKM per Metagenome")
plt.title("Plasmid Abundance per Metagenome by Host Genus (n ≥ 3)")

plt.xticks(rotation=45, ha="right")
plt.tight_layout()

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/boxplot_metagenome_rpkm_by_genus_filtered_final_plasmids.png", dpi=300, bbox_inches="tight")
plt.show()

# Most widespread plasmids

# Compute presence (number of metagenomes with RPKM ≥ 1)
binary_matrix = (filtered_df >= 1).astype(int)
presence_counts = binary_matrix.sum(axis=1)

# Compute total and mean abundance
total_abundance = filtered_df.sum(axis=1)
mean_abundance = filtered_df.replace(0, np.nan).mean(axis=1)

# Define 90th percentiles (top 10%)
presence_thresh = presence_counts.quantile(0.90)
total_abundance_thresh = total_abundance.quantile(0.90)

# Get plasmids above thresholds
most_widespread = presence_counts[presence_counts >= presence_thresh]
top_abundant_total = total_abundance[total_abundance >= total_abundance_thresh]

# Intersect widespread and abundant plasmids
widespread_set = set(most_widespread.index)
abundant_set = set(top_abundant_total.index)
top_both = widespread_set & abundant_set

print(f"{len(top_both)} plasmids are both highly widespread and abundant (top 10%).")

# Save to file
filtered_df.loc[list(top_both)].to_csv(
    "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmids.csv"
)
# Restrict to top abundant and widespread plasmids
top_df = filtered_df.loc[list(top_both)]

# 1. Number of metagenomes with RPKM ≥ 1 per plasmid
top_presence_counts = (top_df >= 1).sum(axis=1)

# 2. Mean RPKM per plasmid (non-zero values only)
top_mean_rpkm = top_df.replace(0, np.nan).mean(axis=1)

# 3. Number of *sponge genera* in which each plasmid is found

# Transpose top_df to get metagenomes as rows
top_df_T = top_df.T

# Merge with genus info
top_df_with_genus = top_df_T.merge(host_info, left_index=True, right_on="Run", how="left")

# Group by plasmid (i.e., column) and count unique genera with RPKM ≥ 1
def count_genera_with_presence(plasmid_column):
    present = top_df_with_genus[top_df_with_genus[plasmid_column] >= 1]
    return present["biome_genus"].nunique()

top_host_genus_counts = pd.Series({plasmid: count_genera_with_presence(plasmid) for plasmid in top_df.columns})

# Combine into summary DataFrame
summary_df = pd.DataFrame({
    "Metagenome_Presence_Count": top_presence_counts,
    "Mean_RPKM": top_mean_rpkm,
    "Host_Genus_Count": top_host_genus_counts
})

# Save summary to CSV
summary_df.to_csv("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_summary.csv")

# Plot
plt.figure(figsize=(10, 7))

# Prepare color and size based on category
colors = ['red' if idx in top_both else 'gray' for idx in clustered_df.index]
sizes = [50 if idx in top_both else 15 for idx in clustered_df.index]

plt.scatter(
    clustered_df["metagenome_presence_count"],
    total_abundance[clustered_df.index] + 1e-3,  # Avoid log(0)
    c=colors,
    s=sizes,
    alpha=0.6,
    edgecolor='k',
    linewidth=0.2
)

# Add threshold lines
plt.axhline(y=total_abundance_thresh, color='darkgray', linestyle='--', label='90th %ile abundance')
plt.axvline(x=presence_thresh, color='darkgray', linestyle='--', label='90th %ile presence')

# Annotate top 5 most abundant among top_both
highlight_abundance = total_abundance.loc[list(top_both)]
top_labels = highlight_abundance.sort_values(ascending=False).head(5)
for plasmid_id in top_labels.index:
    x = clustered_df.loc[plasmid_id, "metagenome_presence_count"]
    y = total_abundance.loc[plasmid_id]
    plt.text(x, y + 0.05 * y, plasmid_id[:20], fontsize=7, ha='center')

# Axis and aesthetics
plt.yscale("log")
plt.xlabel("Number of Metagenomes with RPKM ≥ 1", fontsize=12)
plt.ylabel("Total RPKM across Metagenomes (log scale)", fontsize=12)
plt.title(f"Plasmid Widespreadness vs. Abundance\n({len(top_both)} plasmids in top 10% both)", fontsize=14)
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.3)
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_vs_abundance_highlighted_final_plasmids.png", dpi=300)
plt.show()
# Distribution of top widespread and abundant plasmids across host genera
# Save plasmid IDs to a text file
with open("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_names.txt", "w") as f:
    for plasmid_id in sorted(top_both):
        f.write(plasmid_id + "\n")

# Filter to top_both plasmids
top_both_df = filtered_df.loc[filtered_df.index.intersection(top_both)]

# Count number of metagenomes per genus
genus_map = merged_df[["Metagenome", "biome_genus"]].copy()
genus_map = genus_map.rename(columns={"biome_genus": "Genus"})  # Optional rename for consistency
genus_counts = genus_map["Genus"].value_counts()

# Filter to keep only genera with more than 1 metagenome
valid_genera = genus_counts[genus_counts > 1].index
genus_map = genus_map[genus_map["Genus"].isin(valid_genera)]

# Create genus_with_count column dynamically
genus_map["genus_with_count"] = genus_map["Genus"] + " (n=" + genus_map["Genus"].map(genus_counts).astype(str) + ")"

# Merge it back into merged_df
merged_df = merged_df.drop(columns=["genus_with_count"], errors="ignore")
merged_df = merged_df.merge(genus_map[["Metagenome", "genus_with_count"]], on="Metagenome", how="left")

# Melt to long format and merge with genus info 
long_df = top_both_df.reset_index().melt(
    id_vars="Genome", var_name="Metagenome", value_name="RPKM"
).rename(columns={"Genome": "Plasmid"})

long_df = long_df.merge(merged_df[["Metagenome", "genus_with_count"]], on="Metagenome", how="left")
long_df = long_df[long_df["RPKM"] > 0]

# Count number of metagenomes per genus 
total_per_genus = merged_df[["Metagenome", "genus_with_count"]].drop_duplicates()
metagenomes_per_genus = total_per_genus["genus_with_count"].value_counts()
metagenomes_per_genus = metagenomes_per_genus.reindex(sorted(metagenomes_per_genus.index))

# Create presence matrix and normalize 
# After filtering long_df for RPKM > 0, drop duplicates on Plasmid and Metagenome:
long_df_unique = long_df.drop_duplicates(subset=['Plasmid', 'Metagenome'])

# Then create the presence count matrix
presence_count_matrix = long_df_unique.pivot_table(
    index='Plasmid',
    columns='genus_with_count',
    aggfunc='size',
    fill_value=0
)

metagenomes_per_genus = metagenomes_per_genus.reindex(presence_count_matrix.columns)
relative_presence_matrix = presence_count_matrix.divide(metagenomes_per_genus, axis=1) * 100
print("Summary statistics:")
print(relative_presence_matrix.describe())

# HMA-LMA status dictionary (https://www.nature.com/articles/s41598-018-26641-9, 
#https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.00752/full),
#https://link.springer.com/article/10.1007/s002270000503, 10.3389/fmicb.2021.771589
#https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00877-y
hma_lma_dict = {
    "Coscinoderma": "HMA",
    "Agelas": "HMA",
    "Geodia": "HMA",
    "Amphimedon": "LMA",
    "Cinachyrella": "LMA",
    "Niphates": "LMA",
    "Mycale": "LMA",
    "Spongilla": "LMA",
    "Halichondria": "LMA",
    "Haliclona": "LMA",
    "Cliona": "LMA",
    "Phyllospongia": "LMA",
    "Rhopaloeides": "LMA", 
    "Manihinea": "Unknown",
    "Theonella": "HMA",
    "Aiolochroia": "HMA",
    "Pseudoceratina": "HMA",
    "Xestospongia": "LMA",
    "Aplysina": "HMA",
    "Ircinia": "HMA",
    "Ianthella": "LMA",
    "Isodictya": "LMA",
    "Ephydatia": "LMA",
    "Pericharax": "Unknown",
    "Cymbastela": "LMA",
    "Lamellodysidea": "LMA",
    
}

# Extract genus from genus_with_count
genus_base_names = [label.split(" (n=")[0] for label in relative_presence_matrix.columns]

hma_lma_colors = {
    "HMA": "darkgreen",
    "LMA": "orange",
    "Unknown": "lightgrey"
}

hma_lma_labels = [hma_lma_dict.get(genus, "Unknown") for genus in genus_base_names]
col_colors = [hma_lma_colors.get(label, "lightgrey") for label in hma_lma_labels]


# Plot
relative_presence_matrix.to_csv("C:/Users/hayat/Downloads/R_files/data/top_279_plasmid_relative_presence_matrix_filtered.csv")


g = sns.clustermap(
    relative_presence_matrix,
    cmap="Blues",
    figsize=(14, 10),
    cbar_kws={"label": "Presence (% of metagenomes)"},
    method="average",
    metric="euclidean",
    dendrogram_ratio=(0.1, 0.1),
    colors_ratio=0.01,
    col_colors=col_colors  # Add host HMA/LMA annotation here
)

# Adjust colorbar position and size
pos = g.ax_cbar.get_position()
g.ax_cbar.set_position([pos.x0 + 0.9, pos.y0, pos.width * 0.5, pos.height])

# Continue with label and tick adjustments
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    fontsize=9,
    rotation=90,
    ha="center",
    #weight="light"
)

g.ax_heatmap.tick_params(axis='y', which='both', length=0)
g.ax_heatmap.set_yticklabels([])

plt.suptitle(
    "Clustered Relative Presence of Top 25% Abundant & Widespread Plasmids\nAcross Host Genera (>1 metagenome with plasmid)",
    fontsize=18,
    y=1.05
)
legend_handles = [
    Patch(color=color, label=label) for label, color in hma_lma_colors.items()
]

g.ax_col_dendrogram.legend(
    handles=legend_handles,
    title="Host Microbial Abundance",
    loc="center",
    bbox_to_anchor=(1.2, -3),
    ncol=3
)
g.ax_heatmap.set_xlabel("Host Genus (with metagenome count)", fontsize=12)
g.ax_heatmap.set_ylabel("Plasmid", fontsize=12)

g.savefig("C:/Users/hayat/Downloads/R_files/graphs/top_10_percent_abundant_widespread_relative_presence_clustermap.png", dpi=300)
plt.show()

# Map metagenomes (columns) to their sponge genus
metagenome_to_genus = host_info.set_index("Run")["biome_genus"].to_dict()

# Create a mapping from genus to list of metagenomes
from collections import defaultdict
genus_to_metagenomes = defaultdict(list)
for metagenome, genus in metagenome_to_genus.items():
    if genus:  # Skip NaN or None
        genus_to_metagenomes[genus].append(metagenome)

# Create a new DataFrame: rows = plasmids, columns = sponge genera
plasmid_vs_genus = pd.DataFrame(0, index=plasmid_presence.index, columns=genus_to_metagenomes.keys())

# Fill in the matrix: if a plasmid is present in any metagenome of a genus, mark 1
for genus, metagenomes in genus_to_metagenomes.items():
    valid_metagenomes = [m for m in metagenomes if m in plasmid_presence.columns]
    genus_presence = (plasmid_presence[valid_metagenomes].sum(axis=1) > 0).astype(int)
    plasmid_vs_genus[genus] = genus_presence

# Save to CSV
plasmid_vs_genus.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmids_vs_sponge_genera.csv")


# PKM per HMA/LMA

rpkm_per_metagenome = filtered_df.sum(axis=0).reset_index()
rpkm_per_metagenome.columns = ["Metagenome", "Total_RPKM"]

rpkm_with_genus = rpkm_per_metagenome.merge(host_info, left_on="Metagenome", right_on="Run", how="left")

rpkm_with_genus["Type"] = rpkm_with_genus["biome_genus"].map(hma_lma_dict)

rpkm_with_genus = rpkm_with_genus.dropna(subset=["Type"])

# plt.figure(figsize=(8, 5))
# sns.boxplot(data=rpkm_with_genus, x="Type", y="Total_RPKM", palette="Set2")
# plt.title("Total Plasmid RPKM per Metagenome by Sponge Type (HMA vs LMA)")
# plt.ylabel("Total RPKM")
# plt.xlabel("Sponge Type")
# plt.yscale("log")  # Optional, use if RPKM values vary widely
# plt.tight_layout()
# sns.stripplot(data=rpkm_with_genus, x="Type", y="Total_RPKM", color="black", size=4, jitter=True, alpha=0.5)
# plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/RPKM_per_metagenome_per_hma_lma.png", dpi=300)


# plt.show()

# mann-whitney-u test


hma_data = rpkm_with_genus[rpkm_with_genus["Type"] == "HMA"]["Total_RPKM"]
lma_data = rpkm_with_genus[rpkm_with_genus["Type"] == "LMA"]["Total_RPKM"]
stat, p_value = mannwhitneyu(hma_data, lma_data, alternative="two-sided")

plt.figure(figsize=(8, 5))

# Basic boxplot
sns.boxplot(data=rpkm_with_genus, x="Type", y="Total_RPKM", palette="Set2", showfliers=True)

# Calculate means and medians
grouped = rpkm_with_genus.groupby("Type")["Total_RPKM"]
means = grouped.mean()
medians = grouped.median()

# Overlay mean and median markers
for i, group in enumerate(means.index):
    plt.scatter(i, means[group], color='blue', marker='D', label='Mean' if i == 0 else "", s=60, edgecolor='black', zorder=10)
    plt.scatter(i, medians[group], color='red', marker='o', label='Median' if i == 0 else "", s=60, edgecolor='black', zorder=10)

plt.yscale("log")
plt.title("Total Plasmid RPKM per Metagenome by Sponge Type (HMA vs LMA)")
plt.ylabel("Total RPKM (log scale)")
plt.xlabel("Sponge Type")

# Add legend for mean and median
plt.legend()

x1, x2 = 0, 1
y, h, col = rpkm_with_genus["Total_RPKM"].max() * 1.1, 0.2, 'k'
plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)

if p_value < 0.0001:
    p_text = "****"
elif p_value < 0.001:
    p_text = "***"
elif p_value < 0.01:
    p_text = "**"
elif p_value < 0.05:
    p_text = "*"
else:
    p_text = "ns"

plt.text((x1 + x2) * 0.5, y + h + 0.05, p_text, ha='center', va='bottom', color=col)

plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_RPKM_per_hma_lma.png", dpi=300)

plt.show()

# Calculate plasmid richness per metagenome (number of plasmids with RPKM ≥ 1)
plasmid_richness = (filtered_df >= 1).sum(axis=0).reset_index()
plasmid_richness.columns = ["Metagenome", "Plasmid_Richness"]

# Merge with host genus info
richness_with_genus = plasmid_richness.merge(host_info, left_on="Metagenome", right_on="Run", how="left")

# Assign HMA/LMA status
richness_with_genus["Type"] = richness_with_genus["biome_genus"].map(hma_lma_dict)

# Drop rows with unknown or missing type
richness_with_genus = richness_with_genus.dropna(subset=["Type"])

# Mann-Whitney U test
hma_richness = richness_with_genus[richness_with_genus["Type"] == "HMA"]["Plasmid_Richness"]
lma_richness = richness_with_genus[richness_with_genus["Type"] == "LMA"]["Plasmid_Richness"]
stat_rich, p_rich = mannwhitneyu(hma_richness, lma_richness, alternative="two-sided")


sns.boxplot(data=richness_with_genus, x="Type", y="Plasmid_Richness", palette="Set2", showfliers=True)

# Calculate means and medians
grouped = richness_with_genus.groupby("Type")["Plasmid_Richness"]
means = grouped.mean()
medians = grouped.median()

# Overlay mean and median markers
for i, group in enumerate(means.index):
    plt.scatter(i, means[group], color='blue', marker='D', label='Mean' if i == 0 else "", s=60, edgecolor='black', zorder=10)
    plt.scatter(i, medians[group], color='red', marker='o', label='Median' if i == 0 else "", s=60, edgecolor='black', zorder=10)

plt.title("Plasmid Diversity per Metagenome by Sponge Type (HMA vs LMA)")
plt.ylabel("Number of Distinct Plasmids (RPKM ≥ 1)")
plt.xlabel("Sponge Type")

# Add legend for mean and median
plt.legend()


# Annotate p-value
x1, x2 = 0, 1
y = max(hma_richness.max(), lma_richness.max()) * 1.05
h = 2
plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c='k')

if p_rich < 0.0001:
    p_text = "****"
elif p_rich < 0.001:
    p_text = "***"
elif p_rich < 0.01:
    p_text = "**"
elif p_rich < 0.05:
    p_text = "*"
else:
    p_text = "ns"

plt.text((x1 + x2) * 0.5, y + h + 1, p_text, ha='center', va='bottom', color='k')

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/plasmid_diversity_per_hma_lma.png", dpi=300)
plt.show()

# Summary statistics for plasmid richness by sponge type
richness_summary = richness_with_genus.groupby("Type")["Plasmid_Richness"].describe()
print("\nPlasmid Richness Summary by Sponge Type (HMA vs LMA):")
print(richness_summary)
