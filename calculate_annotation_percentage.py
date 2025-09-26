import pandas as pd

infile = "C:/Users/hayat/Downloads/R_files/data/eggnog_output.emapper.annotations"

# Read file: skip lines starting with '##', but keep the '#query' header
df = pd.read_csv(infile, sep="\t", comment="#", header=None)

# Re-read with proper header: the first non-comment line is the header
with open(infile) as f:
    header = None
    for line in f:
        if line.startswith("#query"):
            header = line.strip().lstrip("#").split("\t")
            break

df = pd.read_csv(infile, sep="\t", comment="#", names=header)

# Rename for convenience
df = df.rename(columns={"#query": "query"})

# Total proteins
total = len(df)

forbidden = {"-", "S", "R"}   # add "R" here if you want to exclude it too

cog_assigned = df["COG_category"].apply(
    lambda x: isinstance(x, str) 
              and x != "-" 
              and not any(c in forbidden for c in x.split(","))
).sum()


# Proteins with KEGG Pathway (exclude "-")
kegg_assigned = df["KEGG_Pathway"].apply(
    lambda x: isinstance(x, str) and x != "-"
).sum()

# Proteins with PFAMs (exclude "-")
pfam_assigned = df["PFAMs"].apply(
    lambda x: isinstance(x, str) and x != "-"
).sum()

# Print results
print(f"Total proteins: {total}")
print(f"Proteins with COG category: {cog_assigned} ({cog_assigned/total*100:.2f}%)")
print(f"Proteins with KEGG pathway: {kegg_assigned} ({kegg_assigned/total*100:.2f}%)")
print(f"Proteins with PFAM: {pfam_assigned} ({pfam_assigned/total*100:.2f}%)")
