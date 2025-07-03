import pandas as pd

rpkm_df = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmids.csv", sep=",")

# Convert all metagenome columns to numeric
for col in rpkm_df.columns[1:]:
    rpkm_df[col] = pd.to_numeric(rpkm_df[col], errors='coerce')

metadata = pd.read_csv("C:/Users/hayat/Downloads/R_files/data/top_widespread_plasmid_metadata.tsv", sep="\t")

rpkm_threshold = 1.0

presence_counts = (rpkm_df.iloc[:, 1:] > 0).sum(axis=1)
rpkm_df['presence_count'] = presence_counts

top_plasmids = rpkm_df.sort_values('presence_count', ascending=False).head(5)

rows = []
for _, plasmid_row in top_plasmids.iterrows():
    plasmid = plasmid_row['Genome']
# Filter only numeric columns for comparison
    numeric_cols = plasmid_row.index.drop(['Genome', 'presence_count'])

    present_runs = [col for col in numeric_cols if plasmid_row[col] > rpkm_threshold]
    for run in present_runs:
        rows.append({'Plasmid': plasmid, 'Run': run})

presence_df = pd.DataFrame(rows)

# Before merging, ensure columns have consistent names and no duplicates
# For example, keep 'Plasmid' and 'Run' in presence_df
presence_df = presence_df[['Plasmid', 'Run']]

# From metadata, keep only columns needed: 'Run', 'Latitude', 'Longitude', 'host_genus', 'geo_cluster' if exists
metadata_subset = metadata[['Run', 'Latitude', 'Longitude', 'host_genus', 'geo_cluster']]

# Merge on 'Run'
merged = presence_df.merge(metadata_subset, on='Run', how='left')

# Check rows with missing metadata
missing_meta = merged[merged['Latitude'].isna() | merged['Longitude'].isna()]
print(f"Rows with missing metadata: {len(missing_meta)}")

# Optionally, drop rows without metadata
merged = merged.dropna(subset=['Latitude', 'Longitude'])

# Check columns now
print(merged.head())


merged = presence_df.merge(metadata, on='Run', how='left')
merged = merged.rename(columns={'Plasmid_x': 'Plasmid'}).drop(columns=['Plasmid_y'])


for plasmid, group_df in merged.groupby('Plasmid'):
    filename = f"C:/Users/hayat/Downloads/R_files/data/{plasmid}_locations.tsv"
    group_df[['Plasmid', 'Run', 'Latitude', 'Longitude', 'host_genus']].to_csv(
        filename,
        sep="\t", index=False
    )
    print(f"Saved {filename} with {len(group_df)} entries.")
