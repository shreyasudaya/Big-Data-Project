import os
from Bio import SeqIO
import numpy as np
from PIL import Image
from collections import Counter

# Step 1: Define input and output directories
fasta_dir = "Dataset/sequences"  # Path to the directory containing individual FASTA files
output_base_dir = "Dataset/kmer3"  # Base directory for output images

# Step 2: Function to extract 3-mers from a sequence
def extract_kmers(seq, k=3):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

# Step 3: Count k-mer frequencies
def kmer_frequencies(seq, k=3):
    kmers = extract_kmers(seq, k)
    kmer_count = Counter(kmers)
    return kmer_count

# Step 4: Convert frequencies to a matrix suitable for image creation
def frequencies_to_image(frequency_dict, img_size=(64, 64)):
    # Create a blank image matrix
    image_matrix = np.zeros((img_size[0], img_size[1]))

    # Normalize frequencies and place them in the image matrix
    max_freq = max(frequency_dict.values(), default=1)  # Avoid division by zero
    for kmer, count in frequency_dict.items():
        index = hash(kmer) % (img_size[0] * img_size[1])  # Use hash for index
        x = index // img_size[1]
        y = index % img_size[1]
        if x < img_size[0] and y < img_size[1]:  # Ensure within bounds
            image_matrix[x, y] = count / max_freq  # Normalize the frequency

    return image_matrix

# Step 5: Process each FASTA file in the specified directory
for fasta_file in os.listdir(fasta_dir):
    if fasta_file.endswith('.fasta'):  # Process only .fasta files
        fasta_path = os.path.join(fasta_dir, fasta_file)
        sequences = SeqIO.parse(fasta_path, "fasta")

        # Process each record in the FASTA file
        for record in sequences:
            # Extract accession number and variant from the header
            header_parts = record.description.split('|')
            accession = header_parts[0].replace('.', '_')  # Replace '.' with '_' in accession for file name
            variant = header_parts[2]

            # Create a directory for the variant if it doesn't exist
            variant_directory = os.path.join(output_base_dir, variant)
            os.makedirs(variant_directory, exist_ok=True)

            # Generate k-mer frequencies and convert to image
            frequency_dict = kmer_frequencies(str(record.seq), k=3)
            img_array = frequencies_to_image(frequency_dict)
            
            # Convert the frequency matrix to an image
            img = Image.fromarray((img_array * 255).astype('uint8'), 'L')  # Create grayscale image

            # Save the image with accession number in the appropriate variant folder
            img.save(os.path.join(variant_directory, f"{accession}.png"))
            print(f"Saved {accession}.png in folder {variant_directory}")

print("All sequences have been processed and images saved.")
