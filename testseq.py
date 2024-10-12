import numpy as np

def create_fcgr_matrix(kmer_counts):
    # Create an empty FCGR matrix of size 64x64
    fcgr_matrix = np.zeros((64, 64))

    # Define a mapping for nucleotide bases to numbers
    base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # Iterate through the kmer_counts to fill the FCGR matrix
    for kmer, count in kmer_counts.items():
        if len(kmer) == 3:  # Only consider 3-mers
            first_base = base_to_num[kmer[0]]
            second_base = base_to_num[kmer[1]]
            third_base = base_to_num[kmer[2]]
            
            # Increment the count in the matrix for the first two bases
            fcgr_matrix[first_base * 4 + second_base][third_base] += count

    return fcgr_matrix

# Example k-mer counts based on your output
kmer_counts = {
    'AAA': 924, 'AAC': 616, 'AAG': 579, 'AAT': 760, 
    'ACA': 809, 'ACC': 375, 'ACG': 165, 'ACT': 674,
    'AGA': 603, 'AGC': 302, 'AGG': 328, 'AGT': 506,
    'ATA': 472, 'ATC': 342, 'ATG': 723, 'ATT': 771,
    'CAA': 702, 'CAC': 458, 'CAG': 437, 'CAT': 486,
    'CCA': 355, 'CCC': 115, 'CCG': 72, 'CCT': 342,
    'CGA': 95, 'CGC': 97, 'CGG': 75, 'CGT': 170,
    'CTA': 556, 'CTC': 287, 'CTG': 493, 'CTT': 741,
    'GAA': 535, 'GAC': 340, 'GAG': 297, 'GAT': 439,
    'GCA': 373, 'GCC': 188, 'GCG': 88, 'GCT': 519,
    'GGA': 280, 'GGC': 223, 'GGG': 134, 'GGT': 454,
    'GTA': 468, 'GTC': 264, 'GTG': 554, 'GTT': 701,
    'TAA': 718, 'TAC': 609, 'TAG': 426, 'TAT': 623,
    'TCA': 546, 'TCC': 206, 'TCG': 112, 'TCT': 542,
    'TGA': 633, 'TGC': 546, 'TGG': 554, 'TGT': 857,
    'TTA': 880, 'TTC': 513, 'TTG': 820, 'TTT': 1007,
}

# Create the FCGR matrix
fcgr_matrix = create_fcgr_matrix(kmer_counts)

# Output the FCGR matrix shape and the matrix itself

print("FCGR matrix shape:", fcgr_matrix.shape)
with open("log.txt", "w") as log_file:
        log_file.write("FCGR matrix shape: {}\n".format(fcgr_matrix.shape))
        log_file.write("FCGR matrix:\n")
        np.savetxt(log_file, fcgr_matrix, fmt='%.5f')
print(fcgr_matrix)
