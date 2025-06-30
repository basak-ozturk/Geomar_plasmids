import pandas as pd
from io import StringIO

# Define the COG categories of interest (metabolic functions)
metabolic_cogs = set("CEFGHIPQ")

input_path = "C:/Users/hayat/Downloads/R_files/data/filtered_foreground_eggnog_annotations.tsv"
output_path = "C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmids_eggnog_metabolic.tsv"

# Step 1: Read the file into memory once and find header line
with open(input_path) as f:
    lines = f.readlines()

# Find header line index and reconstruct a StringIO for pandas
for i, line in enumerate(lines):
    if line.startswith("#query"):
        header_line = i
        break

data_io = StringIO("".join(lines[header_line:]))

# Step 2: Read the data with pandas
df = pd.read_csv(data_io, sep="\t")

# Step 3: Filter for metabolic COG categories
def has_metabolic_cog(cog_str):
    if pd.isna(cog_str):
        return False
    return any(c in metabolic_cogs for c in str(cog_str))

filtered_df = df[df["COG_category"].apply(has_metabolic_cog)]

# Step 4: Save the result
filtered_df.to_csv(output_path, sep="\t", index=False)

print(f"Found {len(filtered_df)} proteins with metabolic COGs.")
