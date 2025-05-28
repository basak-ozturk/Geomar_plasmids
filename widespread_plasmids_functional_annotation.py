import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# Load plasmid names
plasmids = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid_names.csv",
    header=None
)[0].tolist()
plasmid_set = set(plasmids)

# Filter annotations
filtered_lines = []
with open("C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations", "r") as infile, \
     open("C:/Users/hayat/Downloads/R_files/data/filtered_hits_widespread_plasmids.tsv", "w") as outfile:
    
    for line in infile:
        if line.startswith("##"):
            continue
        if line.startswith("#query"):
            outfile.write(line)
            header = line.strip().split('\t')
            continue
        if not line.strip():
            continue

        fields = line.strip().split('\t')
        query_full = fields[0]
        base_id = '_'.join(query_full.split('_')[:-1])
        if base_id in plasmid_set:
            outfile.write(line)
            filtered_lines.append(fields)

df = pd.DataFrame(filtered_lines, columns=header)
df["PFAMs"] = df["PFAMs"].fillna("").astype(str)
df["Description"] = df["Description"].fillna("").astype(str)

# Align PFAMs with their corresponding Descriptions
pfam_entries = []
for _, row in df.iterrows():
    pfams = [pfam.strip() for pfam in row["PFAMs"].split(",") if pfam.strip() and pfam.strip() != "-"]
    for pfam in pfams:
        pfam_entries.append((pfam, row["Description"]))

pfam_desc_df = pd.DataFrame(pfam_entries, columns=["PFAM", "Description"])
pfam_counts = pfam_desc_df["PFAM"].value_counts().reset_index()
pfam_counts.columns = ["PFAM", "Count"]

# Adjust number of top PFAM to be plotted
top_n = 50
top_pfam_df = pfam_counts.head(top_n)
top_pfam_df.to_csv("C:/Users/hayat/Downloads/R_files/data/top_pfam_for_curation.csv", index=False)

top_pfam_df_curated = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_pfam_for_curation_curated_consolidated.csv")

# Optional: Sort by Count for consistent bar order
top_pfam_df_curated = top_pfam_df_curated.sort_values("Count", ascending=False)

# Assign colors by Category
plt.figure(figsize=(14, 10))
sns.set_style("whitegrid")

# Use Seaborn barplot to color by Category
ax = sns.barplot(
    data=top_pfam_df_curated,
    y="PFAM", x="Count",
    hue="Category",
    dodge=False,
    palette="Paired"  
)

# Adjust aesthetics
plt.title("Top PFAM Domains in Widespread Plasmid Proteins (Categorized)", fontsize=16)
plt.xlabel("Count", fontsize=12)
plt.ylabel("PFAM", fontsize=12)
plt.yticks(fontsize=8)
plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_plasmids_top_pfam_plot_categorized.png", dpi=300)
plt.show()

# Count
# Drop missing values, ensure strings, split each string into characters (individual COG codes)
cog_series = df["COG_category"].dropna().astype(str).apply(list).explode()

# Exclude unwanted categories
cog_series = cog_series[~cog_series.isin(["-", "S"])]

# Count occurrences
cog_counts = cog_series.value_counts().reset_index()
cog_counts.columns = ["COG", "Count"]


# Map to functional names
cog_function_map = {
    "J": "Translation",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction",
    "M": "Cell wall/membrane biogenesis",
    "N": "Cell motility",
    #"Z": "Cytoskeleton",
    #"W": "Extracellular structures",
    "U": "Intracellular trafficking",
    "O": "Posttranslational modification",
    "C": "Energy production",
    "G": "Carbohydrate metabolism",
    "E": "Amino acid metabolism",
    "F": "Nucleotide metabolism",
    "H": "Coenzyme metabolism",
    "I": "Lipid metabolism",
    "P": "Inorganic ion transport",
    "Q": "Secondary metabolites biosynthesis",
    "R": "General function prediction",
}

cog_counts["Function"] = cog_counts["COG"].map(cog_function_map)

# Plot
import seaborn as sns
plt.figure(figsize=(10, 8))
sns.barplot(data=cog_counts, y="Function", x="Count", color="skyblue")
plt.xlabel("Number of Proteins")
plt.ylabel("COG Functional Category")
plt.title("COG Categories in Widespread Plasmid Proteins")
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_plasmids_cog_category_plot_filtered.png", dpi=300)
plt.show()

