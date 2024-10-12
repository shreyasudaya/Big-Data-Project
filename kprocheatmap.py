import os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from numba import jit
import seaborn as sns

# JIT-compiled function to generate FCGR for k=3
@jit(nopython=True)
def generate_fcgr(sequence, k=3):
    # Define the grid size for k-mers (4^k = 64 for k=3)
    grid_size = 4 ** k
    fcgr = np.zeros((grid_size, grid_size), dtype=np.float64)

    # Define mapping for nucleotides to numbers
    base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    # Iterate through the sequence to form k-mers and populate FCGR
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        
        # Check if all bases in the k-mer are valid
        is_valid_kmer = True
        for base in kmer:
            if base not in base_to_num:
                is_valid_kmer = False
                break

        if is_valid_kmer:
            x = 0
            y = 0
            for base in kmer:
                x = (x << 1) | (base_to_num[base] >> 1)
                y = (y << 1) | (base_to_num[base] & 1)
            fcgr[x, y] += 1

    # Normalize the FCGR
    total = fcgr.sum()
    if total > 0:
        fcgr /= total

    return fcgr

# Path to directory containing variant folders
base_dir = "Dataset/sequences"
output_base_dir = "Dataset/kmers3"

# Iterate through all subfolders (variants)
for variant_folder in os.listdir(base_dir):
    variant_path = os.path.join(base_dir, variant_folder)
    
    # Ensure the path is a directory
    if os.path.isdir(variant_path):
        # Create output directory for the variant
        output_variant_dir = os.path.join(output_base_dir, variant_folder)
        os.makedirs(output_variant_dir, exist_ok=True)
        
        # Process each fasta file in the variant folder
        fasta_files = [f for f in os.listdir(variant_path) if f.endswith(".fasta")]
        for fasta_file in fasta_files:
            fasta_path = os.path.join(variant_path, fasta_file)
                
            # Read the fasta file using BioPython
            for record in SeqIO.parse(fasta_path, "fasta"):
                header = record.id.split('.')[0]  # Extract accession part
                sequence = str(record.seq).upper()  # Get sequence
                
                # Generate FCGR matrix
                fcgr_matrix = generate_fcgr(sequence, k=3)
                
                # Save FCGR matrix as a PNG image
                output_filename = f"{header}.png"
                output_filepath = os.path.join(output_variant_dir, output_filename)
                
                # Create a heatmap
                plt.figure(figsize=(8, 8))
                sns.heatmap(fcgr_matrix, cmap='viridis', annot=False, fmt=".2f", linewidths=.5)
                plt.title("FCGR Matrix")
                plt.axis('off')  # Hide the axes
                plt.savefig(output_filepath, bbox_inches='tight', pad_inches=0, dpi=300)  # Use higher dpi for better quality
                plt.close()
                print(f"Image created and stored: {output_filepath}")

print("FCGR images generated successfully.")
