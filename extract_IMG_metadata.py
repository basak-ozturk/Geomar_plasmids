# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 09:17:48 2025

@author: hayat
"""
import pandas as pd
#import matplotlib.pyplot as plt

# Load the .m8 file (no header, tab-separated)
m8_cols = [
    "query", "subject", "identity", "aln_length", "mismatches", "gap_opens",
    "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"
]
m8 = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/resultAlignments_GlobalSpongePlasmidome_vs_IMGPR.m8", sep="\t", names=m8_cols)

# Extract plasmid_id from 'subject' (everything before the first "|")
m8["plasmid_id"] = m8["subject"].str.split("|").str[0]

# Load the metadata
meta = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/IMGPR_plasmid_data.tsv", sep="\t")

# Merge on plasmid_id
merged = m8.merge(meta, on="plasmid_id", how="left")

# Select useful columns (e.g., your query plasmid, plasmid_id, source_type, ecosystem)
result = merged[["query", "plasmid_id", "source_type", "ecosystem"]].drop_duplicates()

# Save output
result.to_csv("C:/Users/hayat/Downloads/R_files/data/plasmid_matches_with_metadata.tsv", sep="\t", index=False)
print(result.head())

# Print unique ecosystems and their counts
print("\n=== Ecosystem overview ===")
eco_counts = result["ecosystem"].value_counts(dropna=False)
print(eco_counts)



def clean_ecosystem(value):
    # 1) Replace NaN with "Isolate"
    if pd.isna(value):
        return "Isolate"

    parts = [p.strip() for p in value.split(";")]

    # 2) Host-associated;Porifera;Sponge and Host-associated;Porifera → "Host-associated, Porifera"
    if parts[0] == "Host-associated":
        if len(parts) > 1 and parts[1] == "Porifera":
            return "Host-associated, Porifera"
        else:
            return "Host-associated, other"

    # 3) Engineered → collapse all into "Engineered Environment"
    if parts[0] == "Engineered":
        return "Engineered Environment"

    # 4) Environmental → take the third level if it exists
    if parts[0] == "Environmental":
        if len(parts) >= 3:
            return f"{parts[1]}, {parts[2]}"
        elif len(parts) >= 2:
            return parts[1]
        else:
            return "Environmental (unspecified)"

    # Fallback (if something unexpected shows up)
    return value

# Apply cleaning function
result["ecosystem_clean"] = result["ecosystem"].apply(clean_ecosystem)

# Print summary after cleaning
print("\n=== Cleaned ecosystem overview ===")
print(result["ecosystem_clean"].value_counts())


# Save cleaned ecosystem counts
eco_counts = result["ecosystem_clean"].value_counts().reset_index()
eco_counts.columns = ["ecosystem", "count"]
eco_counts.to_csv("C:/Users/hayat/Downloads/R_files/data/ecosystem_clean_counts.tsv", sep="\t", index=False)

# Save source_type counts
src_counts = result["source_type"].value_counts().reset_index()
src_counts.columns = ["source_type", "count"]
src_counts.to_csv("C:/Users/hayat/Downloads/R_files/data/source_type_counts.tsv", sep="\t", index=False)

# Save host-associated breakdown
host_rows = result[result["ecosystem_clean"].str.startswith("Host-associated")]
host_breakdown = host_rows["ecosystem"].value_counts().reset_index()
host_breakdown.columns = ["ecosystem_original", "count"]
host_breakdown.to_csv("C:/Users/hayat/Downloads/R_files/data/host_associated_breakdown.tsv", sep="\t", index=False)

print("Supplementary tables written:")
print("- ecosystem_clean_counts.tsv")
print("- source_type_counts.tsv")
print("- host_associated_breakdown.tsv")

# Keep only host-associated hits
host_hits = result[result["ecosystem_clean"].str.startswith("Host-associated")]

# Group by query plasmid and check the unique ecosystem_clean categories
query_ecosystems = host_hits.groupby("query")["ecosystem_clean"].unique()

# Select queries that have **only Porifera**
exclusive_porifera = query_ecosystems[query_ecosystems.apply(lambda x: list(x) == ["Host-associated, Porifera"])]

print(f"Number of plasmids exclusively found in Porifera: {len(exclusive_porifera)}")
print("Example plasmids:", exclusive_porifera.index.tolist()[:10])

# Top widespread and abundant plasmids
top_plasmids = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_plasmid_names.txt", 
    header=None, names=["query"]
)

# Filter 'result' to keep only plasmids in your top list
top_matches = result[result["query"].isin(top_plasmids["query"])]

# Check how many top plasmids have IMG hits
num_with_hits = top_matches["query"].nunique()
print(f"{num_with_hits} of the top widespread plasmids have at least one IMG/PR hit.")

# Optional: save these results to a CSV
#top_matches.to_csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmids_with_IMG_hits.tsv", sep="\t", index=False)
