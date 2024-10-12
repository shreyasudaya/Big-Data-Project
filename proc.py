import os
import pandas as pd
from Bio import SeqIO

def one_hot_encode(sequence, fasta_file_path):  # Pass fasta_file_path as argument
    """Encodes a DNA sequence into a one-hot representation.

    Args:
        sequence: A DNA sequence.
        fasta_file_path: Path to the FASTA file containing the sequence.

    Returns:
        A pandas DataFrame representing the one-hot encoded sequence, or None if an error occurs.
    """

    encoding_map = {'A': [1, 0, 0, 0],
                    'C': [0, 1, 0, 0],
                    'G': [0, 0, 1, 0],
                    'T': [0, 0, 0, 1]}

    encoded_sequence = []
    for nucleotide in sequence:
        try:
            encoded_sequence.append(encoding_map[nucleotide])
        except KeyError:
            print(f"Error: Unknown nucleotide '{nucleotide}' encountered in file {fasta_file_path}")
            return None  # Or handle the error differently

    return pd.DataFrame(encoded_sequence, columns=['A', 'C', 'G', 'T'])

def process_fasta_files(sequences_dir, one_hot_dir):
    """Processes FASTA files in the specified directory and saves their one-hot encoded representations.

    Args:
        sequences_dir: The directory containing the FASTA files.
        one_hot_dir: The directory where the one-hot encoded CSV files will be saved.
    """

    if not os.path.exists(one_hot_dir):
        os.makedirs(one_hot_dir)

    for subdir in os.listdir(sequences_dir):
        subdir_path = os.path.join(sequences_dir, subdir)
        if os.path.isdir(subdir_path):
            for fasta_file in os.listdir(subdir_path):
                fasta_file_path = os.path.join(subdir_path, fasta_file)
                if fasta_file.endswith('.fasta'):
                    records = list(SeqIO.parse(fasta_file_path, "fasta"))
                    if len(records) > 0:
                        sequence = records[0].seq.upper()
                        encoded_sequence = one_hot_encode(sequence, fasta_file_path)  # Pass the path
                        if encoded_sequence is not None:
                            csv_file_name = os.path.splitext(fasta_file)[0] + ".csv"
                            csv_file_path = os.path.join(one_hot_dir, subdir, csv_file_name)
                            os.makedirs(os.path.dirname(csv_file_path), exist_ok=True)
                            encoded_sequence.to_csv(csv_file_path, index=False)
                            print(f"Saved {csv_file_name}")

# Example usage
sequences_dir = "Dataset\sequences"
one_hot_dir = "Dataset\one-hot"
process_fasta_files(sequences_dir, one_hot_dir)