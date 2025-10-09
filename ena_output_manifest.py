import pandas as pd
import requests
import time
import os

# --- USER PARAMETERS ---
mapping_file = "C:/Users/hayat/Downloads/R_files/data/plasmid_id_mapping.csv"      
study_accession = "PRJEB98506"            
embl_dir = r"C:\Users\hayat\Downloads\R_files\data\EMBL_FILES"
output_manifest = "C:/Users/hayat/Downloads/R_files/data/webin_analysis_manifest.txt"


# --- LOAD MAPPING ---
mapping = pd.read_csv(mapping_file)  # columns: Original_ID, New_ID
mapping['Run'] = mapping['Original_ID'].apply(lambda x: x.split('_')[0])

# --- FUNCTION TO GET SAMPLE ACCESSION ---
def get_sample_accession(run):
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": run,
        "result": "read_run",
        "fields": "run_accession,sample_accession"
    }
    try:
        r = requests.get(url, params=params, timeout=10)
        lines = r.text.strip().split("\n")
        if len(lines) > 1:
            return lines[1].split("\t")[1]
        else:
            return "NA"
    except Exception as e:
        print(f"Error querying {run}: {e}")
        return "NA"

# --- GET SAMPLE ACCESSIONS ---
unique_runs = mapping['Run'].unique()
run_to_sample = {}
for run in unique_runs:
    run_to_sample[run] = get_sample_accession(run)
    time.sleep(0.2)

mapping['Sample'] = mapping['Run'].map(run_to_sample)

# --- GENERATE MANIFEST ---
with open(output_manifest, 'w') as out:
    for idx, row in mapping.iterrows():
        metagenome = row['Run']
        embl_file = os.path.join(embl_dir, f"{row['New_ID']}.embl")
        block = f"""STUDY       {study_accession}
        
SAMPLE      {row['Sample']}
NAME        {row['New_ID']}
TITLE       Plasmid {row['New_ID']} mined from metagenome {metagenome}
DESCRIPTION Putative plasmid sequence {row['Original_ID']} ({row['New_ID']}) assembled from metagenome {metagenome}
FILE        {embl_file}

"""
        out.write(block)

print(f"Done! Webin-CLI analysis manifest written to {output_manifest}")


# Path to the full manifest generated
full_manifest = r"C:\Users\hayat\Downloads\R_files\data\webin_analysis_manifest.txt"
output_dir = r"C:\Users\hayat\Downloads\R_files\data\manifest_batches"
batch_size = 500  # number of plasmids per batch

os.makedirs(output_dir, exist_ok=True)

# Read the full manifest and split by double newline (each block ends with \n\n)
with open(full_manifest, 'r') as f:
    blocks = f.read().strip().split('\n\n')

# Split into batches
for i in range(0, len(blocks), batch_size):
    batch_blocks = blocks[i:i+batch_size]
    batch_file = os.path.join(output_dir, f"manifest_batch_{i//batch_size + 1}.txt")
    with open(batch_file, 'w') as out:
        out.write('\n\n'.join(batch_blocks) + '\n\n')

print(f"Split {len(blocks)} entries into {len(blocks)//batch_size + 1} batches in {output_dir}")
