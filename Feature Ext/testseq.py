import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

# Function to generate the FCGR matrix
def generate_fcgr(sequence, k):
    dim = 2 ** k
    fcgr_matrix = np.zeros((dim, dim))

    # Convert the sequence to numerical values
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        x, y = 0, 0
        valid_kmer = True
        for j, base in enumerate(kmer):
            if base in mapping:
                bit = mapping[base]
                x = (x << 1) | (bit >> 1)
                y = (y << 1) | (bit & 1)
            else:
                valid_kmer = False
                break
        if valid_kmer:
            fcgr_matrix[x, y] += 1

    # Normalize the matrix
    fcgr_matrix /= np.sum(fcgr_matrix)
    return fcgr_matrix

# Function to save FCGR as an image
def save_fcgr_image(fcgr_matrix, output_path):
    plt.imshow(fcgr_matrix, cmap='hot', interpolation='nearest')
    plt.axis('off')  # Hide axes
    plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
    plt.close()

# Main function to process all FASTA files in subdirectories
def process_fasta_files(input_dir, output_dir, k=8):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Walk through each subdirectory
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fasta") or file.endswith(".fa"):
                file_path = os.path.join(root, file)

                # Create the corresponding subdirectory in the output
                relative_path = os.path.relpath(root, input_dir)
                output_subdir = os.path.join(output_dir, relative_path)
                if not os.path.exists(output_subdir):
                    os.makedirs(output_subdir)

                # Read the sequence from the FASTA file
                fasta_sequences = SeqIO.parse(open(file_path), 'fasta')
                for fasta in fasta_sequences:
                    sequence = str(fasta.seq)

                    # Generate the FCGR matrix
                    fcgr_matrix = generate_fcgr(sequence.upper(), k)

                    # Save the FCGR matrix as an image
                    accession = fasta.id
                    output_path = os.path.join(output_subdir, f"{accession}.png")
                    save_fcgr_image(fcgr_matrix, output_path)
                    print(f"Saved FCGR image for {file} in {output_path}")

# Example usage
input_directory = "Dataset\sequences"
output_directory = "Dataset\Fcgr"
process_fasta_files(input_directory, output_directory, k=8)
