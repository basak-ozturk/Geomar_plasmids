import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Step 1: Load plasmid names
plasmids = pd.read_csv(
    "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid_names.csv",
    header=None
)[0].tolist()
plasmid_set = set(plasmids)

# Step 2: Filter annotations
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

# Step 3: Create DataFrame
df = pd.DataFrame(filtered_lines, columns=header)
df["PFAMs"] = df["PFAMs"].fillna("").astype(str)
df["Description"] = df["Description"].fillna("").astype(str)

# Step 4: Align PFAMs with their corresponding Descriptions
pfam_entries = []
for _, row in df.iterrows():
    pfams = [pfam.strip() for pfam in row["PFAMs"].split(",") if pfam.strip() and pfam.strip() != "-"]
    for pfam in pfams:
        pfam_entries.append((pfam, row["Description"]))

# Step 5: Convert to DataFrame and count
pfam_desc_df = pd.DataFrame(pfam_entries, columns=["PFAM", "Description"])
pfam_counts = pfam_desc_df["PFAM"].value_counts().reset_index()
pfam_counts.columns = ["PFAM", "Count"]

# Step 6: Save and Plot
top_n = 50
top_pfam_df = pfam_counts.head(top_n)
top_pfam_df.to_csv("C:/Users/hayat/Downloads/R_files/data/top_pfam_for_curation.csv", index=False)

plt.figure(figsize=(12, 10))
top_pfam_df.plot(kind='barh', x='PFAM', y='Count', legend=False, color='skyblue')
plt.gca().invert_yaxis()
plt.xlabel("Number")
plt.title(f"Top {top_n} PFAM Domains in Widespread Plasmid Proteins")
plt.xticks(fontsize=5)
plt.yticks(fontsize=5)
plt.tight_layout()
plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/widespread_plasmids_top_pfam_plot.png", dpi=300)
plt.show()
