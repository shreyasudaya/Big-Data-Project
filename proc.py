import os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed

# Function to generate one-hot encoding for a nucleotide sequence
def one_hot_encode(sequence):
    base_to_row = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    one_hot_matrix = np.zeros((4, len(sequence)), dtype=np.float32)
    indices = [base_to_row[base] for base in sequence if base in base_to_row]
    one_hot_matrix[indices, np.arange(len(indices))] = 1
    return one_hot_matrix

# Function to process a single FASTA file
def process_fasta(fasta_path, output_variant_dir):
    # Read the fasta file using BioPython
    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.id.split('.')[0]
        sequence = str(record.seq).upper()

        # Generate one-hot encoded matrix
        one_hot_matrix = one_hot_encode(sequence)
        N = one_hot_matrix.shape[1]

        # Save one-hot matrix as a PNG image
        output_filename = f"{header}.png"
        output_filepath = os.path.join(output_variant_dir, output_filename)

        # Plot and save the one-hot encoded matrix as an image
        fig, ax = plt.subplots(figsize=(N / 100, 4))
        ax.imshow(one_hot_matrix, cmap='gray', aspect='auto')
        ax.axis('off')

        # Save the image
        plt.savefig(output_filepath, bbox_inches='tight', pad_inches=0)
        plt.close()

        return output_filepath  # Return the path of the saved image

# Path to directory containing variant folders
base_dir = "Dataset/sequences"
output_base_dir = "Dataset/one-hot"

# Create a ThreadPoolExecutor to process files in parallel
with ThreadPoolExecutor() as executor:
    futures = []

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

                # Submit the processing function to the executor
                future = executor.submit(process_fasta, fasta_path, output_variant_dir)
                futures.append(future)

    # Collect results and print completion messages
    for future in as_completed(futures):
        output_filepath = future.result()
        print(f"Image created and stored: {output_filepath}")

print("One-hot encoded images of size 4xN generated successfully.")
