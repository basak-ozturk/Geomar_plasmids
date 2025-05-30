import pandas as pd
from collections import Counter

# File paths
summary_file = r"C:\Users\hayat\Downloads\R_files\data\hmmer_results\summary_hits.tsv"
conjugative_output_file = r"C:\Users\hayat\Downloads\R_files\data\hmmer_results\conjugative_plasmids.tsv"
original_list_file = r"C:\Users\hayat\Downloads\R_files\data\unique_CONJscan_hit_ids.txt"
in_both_output_file = r"C:\Users\hayat\Downloads\R_files\data\conjugative_in_both.txt"

# === Load and preprocess summary data ===
df = pd.read_csv(summary_file, sep="\t")
# Extract plasmid name (everything before last underscore)
df["plasmid"] = df["hit_id"].apply(lambda x: "_".join(x.split("_")[:-1]))

# === Define relaxase and T4SS flags ===
relaxase_families = ["MOBF", "MOBP", "MOBQ", "MOBH", "MOBT", "MOBC", "MOBV"]
relaxase_pattern = "|".join([f"T4SS_{mob}" for mob in relaxase_families])

df["has_relaxase"] = df["gene_name"].str.contains(relaxase_pattern)
df["has_t4ss"] = df["gene_name"].str.startswith("T4SS_")

# === Group by plasmid to find conjugative plasmids ===
grouped = df.groupby("plasmid").agg(
    relaxase_found=("has_relaxase", "any"),
    t4ss_found=("has_t4ss", "any")
)

conjugative_plasmids_df = grouped[(grouped["relaxase_found"]) & (grouped["t4ss_found"])]

# Save conjugative plasmids to file
conjugative_plasmids_df.to_csv(conjugative_output_file, sep="\t")
print(f"{len(conjugative_plasmids_df)} conjugative plasmids found. Saved to {conjugative_output_file}")

# === Summary stats ===
total_plasmids = df["plasmid"].nunique()
filtered_plasmids = conjugative_plasmids_df.shape[0]
print(f"Total plasmids: {total_plasmids}")
print(f"Conjugative plasmids (relaxase + T4SS): {filtered_plasmids}")
print(f"Eliminated plasmids: {total_plasmids - filtered_plasmids}")

# === Plasmid counts ===
plasmid_counts = Counter(df["plasmid"])
unique_plasmids_reconstructed = len(plasmid_counts)
plasmids_with_multiple_hits = [p for p, c in plasmid_counts.items() if c > 1]

print(f"Total plasmid entries in reconstructed summary: {df.shape[0]}")
print(f"Unique plasmid names in reconstructed summary: {unique_plasmids_reconstructed}")
print(f"Plasmid names appearing more than once: {len(plasmids_with_multiple_hits)}")

# === Load original plasmid list ===
with open(original_list_file) as f:
    original_plasmids = set(line.strip() for line in f if line.strip())

reconstructed_plasmids = set(df["plasmid"])

# Compare original vs reconstructed plasmid sets
new_plasmids = reconstructed_plasmids - original_plasmids
missing_plasmids = original_plasmids - reconstructed_plasmids

print(f"Plasmids in reconstructed but not original: {len(new_plasmids)}")
print(f"Plasmids in original but not reconstructed: {len(missing_plasmids)}")

print("Sample original plasmid names:", list(original_plasmids)[:5])
print("Sample reconstructed plasmid names:", list(reconstructed_plasmids)[:5])

# === Compare conjugative plasmids to original list ===
conjugative_df = pd.read_csv(conjugative_output_file, sep="\t", index_col=0)
conjugative_plasmids = set(conjugative_df.index)

in_both = conjugative_plasmids & original_plasmids
not_in_original = conjugative_plasmids - original_plasmids

print(f"\nConjugative plasmids: {len(conjugative_plasmids)}")
print(f"Found in original list: {len(in_both)}")
print(f"Not found in original list: {len(not_in_original)}")

if not_in_original:
    print("Example plasmids not in original list:")
    for p in list(not_in_original)[:10]:
        print(" ", p)
else:
    print("All conjugative plasmids are present in the original list.")

# === Save plasmids found in both to a text file ===
with open(in_both_output_file, "w") as f_out:
    for plasmid in sorted(in_both):
        f_out.write(plasmid + "\n")

print(f"Plasmids found in both lists saved to {in_both_output_file}")

# === Find conjugative plasmids with integrases ===


eggnog_file = r"C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations"

# Read file, keep only header (#) and data lines, drop ## comment lines
with open(eggnog_file, "r") as f:
    lines = f.readlines()

# Filter out lines starting with '##' (double hash)
filtered_lines = [line for line in lines if not line.startswith("##")]

# Remove leading '#' from the header line (assumed to be the first line starting with '#')
for i, line in enumerate(filtered_lines):
    if line.startswith("#"):
        filtered_lines[i] = line.lstrip("#").strip() + "\n"
        break

# Create a single string to load into pandas
data_str = "".join(filtered_lines)

from io import StringIO
df_egg = pd.read_csv(StringIO(data_str), sep="\t")

# Now 'query' column should be there
def strip_last_underscore(name):
    parts = name.rsplit("_", 1)
    return parts[0] if len(parts) == 2 else name

df_egg['plasmid_base'] = df_egg['query'].apply(strip_last_underscore)

# Assuming you have your 'in_both' set from before
filtered = df_egg[
    (df_egg['plasmid_base'].isin(in_both)) &
    (df_egg['Description'].str.contains("integrase", case=False, na=False))
]

print(f"Found {len(filtered)} rows matching conjugative plasmids with integrase annotation.")
print(filtered[['query', 'Description', 'Preferred_name']].head(10))

output_integrase_file = r"C:/Users/hayat/Downloads/R_files/data/conjugative_with_integrase.tsv"
filtered.to_csv(output_integrase_file, sep="\t", index=False)
print(f"Saved filtered data to {output_integrase_file}")

# Extract unique plasmid base names from filtered rows
unique_plasmids = filtered['query'].apply(strip_last_underscore).unique()

# Save to a text file, one plasmid per line
unique_plasmids_file = r"C:/Users/hayat/Downloads/R_files/data/conjugative_plasmids_with_integrase.txt"
with open(unique_plasmids_file, "w") as f:
    for plasmid in unique_plasmids:
        f.write(plasmid + "\n")

print(f"Saved {len(unique_plasmids)} unique plasmids with integrase to {unique_plasmids_file}")



# Load unique plasmids with integrase
unique_plasmids_file = r"C:/Users/hayat/Downloads/R_files/data/conjugative_plasmids_with_integrase.txt"
with open(unique_plasmids_file) as f:
    unique_plasmids = set(line.strip() for line in f if line.strip())

# Load top abundant plasmid names
top_abundant_file = r"C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid_names.csv"
top_abundant_df = pd.read_csv(top_abundant_file)

# Convert to a set (ensuring it's strings and whitespace-trimmed)
top_plasmid_names = set(top_abundant_df.iloc[:, 0].astype(str).str.strip())

# Find intersection
common_plasmids = unique_plasmids & top_plasmid_names

print(f"Number of plasmids in both lists: {len(common_plasmids)}")
print("Examples:", list(common_plasmids)[:17])

# Optional: save to file
output_file = r"C:/Users/hayat/Downloads/R_files/data/common_integrase_abundant_plasmids.txt"
with open(output_file, "w") as f:
    for plasmid in sorted(common_plasmids):
        f.write(plasmid + "\n")

print(f"Saved common plasmid list to {output_file}")