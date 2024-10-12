import os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

# Function to generate FCGR for k=3
def generate_fcgr(sequence, k=3):
    # Define the grid size for k-mers (4^k = 64 for k=3)
    grid_size = 2 ** k
    fcgr = np.zeros((grid_size, grid_size))
    
    # Define mapping for nucleotides to numbers
    base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    # Iterate through the sequence to form k-mers and populate FCGR
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if all(base in base_to_num for base in kmer):  # Ensure valid bases
            x = y = 0
            for j, base in enumerate(kmer):
                x = (x << 1) | (base_to_num[base] >> 1)
                y = (y << 1) | (base_to_num[base] & 1)
            fcgr[x, y] += 1

    # Normalize the FCGR
    if fcgr.sum() > 0:
        fcgr /= fcgr.sum()

    return fcgr

# Path to directory containing variant folders
base_dir = "Dataset\sequences"
output_base_dir = "Dataset\kmers3"

# Iterate through all subfolders (variants)
for variant_folder in os.listdir(base_dir):
    variant_path = os.path.join(base_dir, variant_folder)
    
    # Ensure the path is a directory
    if os.path.isdir(variant_path):
        # Create output directory for the variant
        output_variant_dir = os.path.join(output_base_dir, variant_folder)
        os.makedirs(output_variant_dir, exist_ok=True)
        
        # Process each fasta file in the variant folder
        for fasta_file in os.listdir(variant_path):
            if fasta_file.endswith(".fasta"):
                fasta_path = os.path.join(variant_path, fasta_file)
                
                # Read the fasta file using BioPython
                for record in SeqIO.parse(fasta_path, "fasta"):
                    header = record.id.split('.')[0]  # Extract accession part
                    sequence = str(record.seq).upper()  # Get sequence
                    
                    # Generate FCGR matrix
                    fcgr_matrix = generate_fcgr(sequence, k=3)
                    
                    # Normalize FCGR matrix (0 to 1)
                    fcgr_matrix = fcgr_matrix / np.max(fcgr_matrix)
                    
                    # Save FCGR matrix as a PNG image
                    output_filename = f"{header}.png"
                    output_filepath = os.path.join(output_variant_dir, output_filename)
                    
                    # Plot and save the FCGR as an image
                    plt.imshow(fcgr_matrix, cmap='gray', interpolation='nearest')
                    #
                    plt.axis('off')
                    plt.savefig(output_filepath, bbox_inches='tight', pad_inches=0)
                    plt.close()
                    print(f"Image created and stored: {output_filepath}")

print("FCGR images generated successfully.")
