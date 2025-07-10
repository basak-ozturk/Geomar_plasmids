# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 15:11:41 2025

@author: hayat
"""

input_file = "C:/Users/hayat/Downloads/R_files/data/itol_integrase.txt"
output_file = "C:/Users/hayat/Downloads/R_files/data/itol_integrase_symbol.txt"

# Read the input file
with open(input_file, "r") as f:
    lines = f.readlines()

symbol_header = [
    "DATASET_SYMBOL\n",
    "SEPARATOR TAB\n",
    "DATASET_LABEL\tIntegrase\n",
    "COLOR\t#000000\n",
    "LEGEND_TITLE\tIntegrase\n",
    "LEGEND_SHAPES\t1\t1\n",  # same shape: circle
    "LEGEND_COLORS\t#000000\t#000000\n",
    "LEGEND_LABELS\tIntegrase+\tIntegrase-\n",
    "DATA\n"
]

symbol_data = []
data_section = False
for line in lines:
    if line.strip() == "DATA":
        data_section = True
        continue
    if not data_section or line.strip() == "":
        continue

    parts = line.strip().split('\t')
    if len(parts) != 3:
        continue

    plasmid_id, _, label = parts
    fill = "1" if label == "Integrase+" else "0"  # 1 = filled, 0 = empty
    shape = "3"  
    size = "20"
    color = "#000000"
    position = "-1"
    symbol_line = f"{plasmid_id}\t{shape}\t{size}\t{color}\t{fill}\t{position}\t{label}\n"
    symbol_data.append(symbol_line)

with open(output_file, "w") as f:
    f.writelines(symbol_header + symbol_data)

print(f"iTOL external symbol file written to: {output_file}")