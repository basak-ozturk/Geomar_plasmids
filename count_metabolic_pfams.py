import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# === File Paths ===
all_path = "C:/Users/hayat/Downloads/R_files/data/all_plasmids_eggnog_metabolic.tsv"
top_path = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmids_eggnog_metabolic.tsv"
output_dir = Path("C:/Users/hayat/Downloads/R_files/graphs")
output_dir.mkdir(parents=True, exist_ok=True)

# === Load Data ===
df_all = pd.read_csv(all_path, sep="\t")
df_top = pd.read_csv(top_path, sep="\t")

# === PFAM Counts (excluding "-") ===
# === PFAM Counts (excluding "-" and unwanted PFAMs) ===
def get_pfam_counts(df):
    unwanted = {"-", "TraG_N", "TraF"}
    pfams = df["PFAMs"].dropna().str.split(",").explode().str.strip()
    pfams = pfams[~pfams.isin(unwanted)]
    return pfams.value_counts()


pfam_counts_all = get_pfam_counts(df_all)
pfam_counts_top = get_pfam_counts(df_top)

# === Plot Top 10 PFAMs ===
def plot_top_pfams(pfam_counts, title, filename):
    top10 = pfam_counts.head(10)
    plt.figure(figsize=(8, 5))
    sns.barplot(x=top10.values, y=top10.index, palette="viridis")
    plt.title(title)
    plt.xlabel("Count")
    plt.ylabel("Pfam")
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.show()
    plt.close()

plot_top_pfams(pfam_counts_all, "Top 10 Metabolic PFAMs in All Plasmids", "top10_pfams_all.png")
plot_top_pfams(pfam_counts_top, "Top 10 Metabolic PFAMs in Top Widespread Plasmids", "top10_pfams_top.png")

# === COG Category Counts (filtered to only metabolic categories) ===
metabolic_cogs = {
    "C": "Energy",
    "G": "Carbohydrates",
    "E": "Amino acids",
    "F": "Nucleotides",
    "H": "Coenzymes",
    "I": "Lipids",
    "P": "Inorganic ions",
    "Q": "Secondary metabolites"
}

def get_filtered_cog_counts(df):
    cogs = df["COG_category"].dropna().astype(str).str.upper().str.replace(" ", "")
    cogs = cogs.str.extractall(r"([A-Z])")[0]
    cogs = cogs[cogs.isin(metabolic_cogs.keys())]
    return cogs.value_counts().sort_values(ascending=False)

def plot_cogs(cog_counts, title, filename):
    # Relabel to category names
    cog_counts.index = [f"[{k}] {metabolic_cogs[k]}" for k in cog_counts.index]
    plt.figure(figsize=(10, 5))
    sns.barplot(x=cog_counts.values, y=cog_counts.index, palette="mako")
    plt.title(title)
    plt.xlabel("Count")
    plt.ylabel("COG Category")
    plt.tight_layout()
    plt.savefig(output_dir / filename, dpi=300)
    plt.show()
    plt.close()

cog_counts_all = get_filtered_cog_counts(df_all)
cog_counts_top = get_filtered_cog_counts(df_top)

plot_cogs(cog_counts_all, "Metabolic COGs in All Plasmids", "metabolic_cogs_all.png")
plot_cogs(cog_counts_top, "Metabolic COGs in Top Widespread Plasmids", "metabolic_cogs_top.png")

print("Done. Filtered metabolic COG plots saved in:", output_dir)
