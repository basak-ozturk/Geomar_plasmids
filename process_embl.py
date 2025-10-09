import os
import gzip

# --- USER PARAMETERS ---
input_dir = r"C:\Users\hayat\Downloads\R_files\data\EMBL_FILES"
output_dir = r"C:\Users\hayat\Downloads\R_files\data\EMBL_FIXED"
organism_name = "uncultured bacterium"
organism_classification = "Bacteria"

os.makedirs(output_dir, exist_ok=True)

def fix_embl_gz(in_path, out_path):
    with gzip.open(in_path, 'rt') as f:
        lines = f.readlines()

    fixed_lines = []
    in_source = False
    source_fixed = False

    for i, line in enumerate(lines):
        # --- Fix LOCUS line ---
        if line.startswith("ID"):
            parts = line.strip().split(";")
            while len(parts) < 6:  # ensure enough fields
                parts.append("")
            parts[2] = " linear"
            parts[3] = " genomic DNA"
            line = ";".join(parts) + "\n"
            fixed_lines.append(line)
            continue

        # --- Fix OS/OC header lines ---
        if line.startswith("OS   ."):
            fixed_lines.append(f"OS   {organism_name}\n")
            continue
        if line.startswith("OC   ."):
            fixed_lines.append(f"OC   {organism_classification}.\n")
            continue

        # --- Remove /locus_tag lines ---
        if "/locus_tag=" in line:
            continue

        # --- Handle source feature ---
        if line.startswith("FT   source"):
            in_source = True
            source_fixed = False
            fixed_lines.append(line)
            continue

        if in_source:
            if not source_fixed:
                # Insert organism and mol_type immediately after FT source
                fixed_lines.append(f'FT                   /organism="{organism_name}"\n')
                fixed_lines.append('FT                   /mol_type="genomic DNA"\n')
                source_fixed = True
            # End of source block if next feature starts
            if line.startswith("FT   ") and not line.startswith("FT   source"):
                in_source = False

        # --- Skip stray /mol_type outside source ---
        if "/mol_type=" in line and not line.startswith("FT   source"):
            continue

        fixed_lines.append(line)

    # --- Write fixed file as .gz ---
    with gzip.open(out_path, 'wt') as out:
        out.writelines(fixed_lines)

# --- Process all .embl.gz files ---
embl_gz_files = [f for f in os.listdir(input_dir) if f.lower().endswith(".embl.gz")]

if not embl_gz_files:
    print("No .embl.gz files found in input folder!")
else:
    print(f"Found {len(embl_gz_files)} files, processing...")

for filename in embl_gz_files:
    in_file = os.path.join(input_dir, filename)
    out_file = os.path.join(output_dir, filename)
    fix_embl_gz(in_file, out_file)
    print(f"Processed {filename} -> {out_file}")

print("All files processed. Files should now validate in ENA Webin.")
