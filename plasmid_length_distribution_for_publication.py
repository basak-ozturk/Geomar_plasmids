# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 09:39:56 2025

@author: hayat
"""

import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

def analyze_fasta(fasta_file):
    # Read sequences and calculate lengths
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    
    # Calculate statistics
    min_length = min(lengths)
    max_length = max(lengths)
    median_length = np.median(lengths)
    
    # Print statistics
    print(f"Minimum length: {min_length}")
    print(f"Maximum length: {max_length}")
    print(f"Median length: {median_length}")
    
    # Plot histogram
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=50, edgecolor='black')
    plt.title("Distribution of Sequence Lengths for All Putative Sponge Plasmids")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    
    # Add vertical lines for min, max, and median
    plt.axvline(min_length, color='r', linestyle='dashed', linewidth=2, label=f'Min: {min_length}')
    plt.axvline(max_length, color='g', linestyle='dashed', linewidth=2, label=f'Max: {max_length}')
    plt.axvline(median_length, color='b', linestyle='dashed', linewidth=2, label=f'Median: {median_length}')
    
    plt.legend()
    plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/sequence_length_distribution_all_sponge_plas.png", dpi=300)
    plt.close()

analyze_fasta('C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.fasta')
