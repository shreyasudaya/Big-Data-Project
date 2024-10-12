import os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

# Function to generate FCGR for k=3
def generate_fcgr(sequence, k=3):
    # Define the grid size for k-mers (4^k = 64 for k=3)
    grid_size = 4 ** k
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

# Path to the FASTA file
fasta_path = "testsequence.fasta"  # Assuming it's in the same directory

# Read the fasta file using BioPython
for record in SeqIO.parse(fasta_path, "fasta"):
    header = record.id.split('.')[0]  # Extract accession part
    sequence = str(record.seq).upper()  # Get sequence
    
    # Generate FCGR matrix
    fcgr_matrix = generate_fcgr(sequence, k=3)
    
    # Normalize FCGR matrix (0 to 1)
    if np.max(fcgr_matrix) > 0:
        fcgr_matrix = fcgr_matrix / np.max(fcgr_matrix)

    # Save FCGR matrix as a PNG image
    output_filepath = "testkprocoutput.png"  # Specify output filename
    print(fcgr_matrix)
    # Plot and save the FCGR as an image
    plt.imshow(fcgr_matrix, cmap='gray', interpolation='nearest')
    plt.axis('off')
    plt.savefig(output_filepath, bbox_inches='tight', pad_inches=0)
    plt.close()
    print(f"Image created and stored: {output_filepath}")

print("FCGR image generated successfully.")
