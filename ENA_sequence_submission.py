# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 08:24:07 2025

@author: hayat
"""

import os
import subprocess
#from datetime import datetime

# === CONFIGURATION ===

manifest_dir = r"C:\Users\hayat\Downloads\R_files\data\manifest_fixed"
webin_cli_path = r"C:\Users\hayat\Downloads\webin-cli-9.0.1.jar"
username = "boezturk@geomar.de"
password = ""
mode = "submit" 
log_dir = r"C:\Users\hayat\Downloads\R_files\data\webin_logs"
# === SETUP ===
os.makedirs(log_dir, exist_ok=True)

# Get manifest files
manifest_files = sorted([f for f in os.listdir(manifest_dir) if f.endswith(".txt")])
if not manifest_files:
    raise FileNotFoundError("No manifest files found in manifest_fixed/")

# === MAIN LOOP ===
for manifest in manifest_files:
    manifest_path = os.path.join(manifest_dir, manifest)
    base_name = os.path.splitext(manifest)[0]
    this_log_dir = os.path.join(log_dir, base_name)
    os.makedirs(this_log_dir, exist_ok=True)
    
    print(f"\n=== Submitting {manifest} ===")
    
    report_file = os.path.join(this_log_dir, "webin-cli.report")
    
    # Skip if already submitted
    if os.path.exists(report_file) and "analysis accession" in open(report_file).read().lower():
        print(f"→ Skipping {manifest}, already submitted.")
        continue

    # Build command
    cmd = [
        "java", "-jar", webin_cli_path,
        #"-mode", mode,
        "-context", "sequence",
        "-manifest", manifest_path,
        "-outputDir", this_log_dir,
        "-username", username,
        "-password", password
    ]
    if mode == "submit":
        cmd.append("-submit")
    
    # Run submission
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Write stdout/stderr logs
    with open(os.path.join(this_log_dir, "webin_output.log"), "w") as log_out:
        log_out.write(result.stdout)
    with open(os.path.join(this_log_dir, "webin_error.log"), "w") as log_err:
        log_err.write(result.stderr)
    
    print(result.stdout)
    if result.returncode == 0:
        print(f"✅ Successfully processed {manifest}")
    else:
        print(f"❌ Error in {manifest} — check {this_log_dir}/webin_error.log")

print("\nAll manifests processed.")
