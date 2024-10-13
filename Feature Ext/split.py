import os
from Bio import SeqIO

# Step 1: Define input and output directories
fasta_file = "sequences.fasta"  # Path to your input FASTA file
output_base_dir = "Dataset\sequences"  # Base directory for output

# Step 2: Create the output base directory if it doesn't exist
if not os.path.exists(output_base_dir):
    os.makedirs(output_base_dir)

# Step 3: Read the FASTA file and process each sequence
for record in SeqIO.parse(fasta_file, "fasta"):
    # Extract accession number and variant from the header
    header_parts = record.description.split('|')
    accession = header_parts[0].replace('.', '_')  # Replace '.' with '_' in accession for file name
    variant = header_parts[2]

    # Create a directory for the variant if it doesn't exist
    variant_dir = os.path.join(output_base_dir, variant)
    if not os.path.exists(variant_dir):
        os.makedirs(variant_dir)

    # Define the output file path for the individual FASTA file
    output_file_path = os.path.join(variant_dir, f"{accession}.fasta")

    # Write the sequence to the individual FASTA file
    with open(output_file_path, "w") as output_file:
        SeqIO.write(record, output_file, "fasta")

    print(f"Saved {accession}.fasta in folder {variant_dir}")

print("All sequences have been processed and saved.")
