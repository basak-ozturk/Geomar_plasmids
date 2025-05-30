import os
import glob

# Path to the folder with .res_hmm_extract files
folder_path = r"C:\Users\hayat\Downloads\R_files\data\hmmer_results"

# Output file path
output_file = "C:/Users/hayat/Downloads/R_files/data/summary_hits.tsv"

# Expected column names
header = [
    "hit_id", "replicon_name", "position_hit", "hit_sequence_length", 
    "gene_name", "i_eval", "score", "profile_coverage", 
    "sequence_coverage", "begin", "end"
]

summary_lines = []

# Loop over all extract files
for file_path in glob.glob(os.path.join(folder_path, "*.res_hmm_extract")):
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # skip comments and empty lines

            parts = line.split("\t")
            if len(parts) < 11:
                print(f"⚠️ Skipping malformed line in {os.path.basename(file_path)}:\n{line}")
                continue

            summary_lines.append(line)

# Write output
if summary_lines:
    with open(output_file, "w", encoding="utf-8") as out_f:
        out_f.write("\t".join(header) + "\n")
        for line in summary_lines:
            out_f.write(line + "\n")
    print(f"✅ Summary written to: {output_file}")
else:
    print("❌ No valid lines found. Nothing written.")