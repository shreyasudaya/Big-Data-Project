import os

def get_file_stems(directory, extension):
    """
    Get a set of file stems (filenames without extension) for all files in the given directory
    and its subdirectories, filtering by the given extension.
    """
    file_stems = {}
    
    for variant in os.listdir(directory):
        variant_path = os.path.join(directory, variant)
        if os.path.isdir(variant_path):
            file_stems[variant] = set(
                os.path.splitext(f)[0] for f in os.listdir(variant_path) if f.endswith(extension)
            )
    
    return file_stems

def delete_unmatched_files(sequences_dir, onehot_dir):
    """
    Deletes FASTA files in `sequences_dir` that do not have corresponding CSV files in `onehot_dir`.
    """
    # Get file stems for both directories
    fasta_files = get_file_stems(sequences_dir, '.fasta')
    csv_files = get_file_stems(onehot_dir, '.csv')
    
    # Loop through variants and delete unmatched FASTA files
    for variant, fasta_stems in fasta_files.items():
        if variant in csv_files:
            unmatched_files = fasta_stems - csv_files[variant]  # Files in sequences but not in one-hot
            
            variant_path = os.path.join(sequences_dir, variant)
            for file_stem in unmatched_files:
                fasta_file = os.path.join(variant_path, file_stem + '.fasta')
                os.remove(fasta_file)
                print(f"Deleted: {fasta_file}")
        else:
            print(f"Warning: No matching subfolder found for variant {variant} in one-hot directory")

# Paths
sequences_directory = "Dataset\sequences"  # Replace with the actual path
onehot_directory = "Dataset\one-hot"       # Replace with the actual path

# Delete unmatched files
delete_unmatched_files(sequences_directory, onehot_directory)
