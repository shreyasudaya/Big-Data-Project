from pyspark.sql import SparkSession
from pyspark.sql.functions import udf
from pyspark.sql.types import ArrayType, IntegerType
from Bio import SeqIO
import numpy as np
from PIL import Image
import os


# Initialize Spark session
spark = SparkSession.builder \
    .appName("FASTA Sequence Encoding") \
    .getOrCreate()

# Function to read FASTA files and return a DataFrame
def read_fasta(fasta_file):
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        records.append((record.id, str(record.seq), record.description))
    return spark.createDataFrame(records, ["accession", "sequence", "description"])

# One-hot encoding function
def one_hot_encode(seq):
    mapping = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'C': [0, 0, 0, 1]}
    return [mapping[nuc] for nuc in seq if nuc in mapping]

# Register the UDF for one-hot encoding
one_hot_udf = udf(one_hot_encode, ArrayType(ArrayType(IntegerType())))

# Path to the directory containing the FASTA files
fasta_dir = "Dataset/sequences"  # Modify this path as needed

# List all FASTA files in the directory
fasta_files = [os.path.join(fasta_dir, file) for file in os.listdir(fasta_dir) if file.endswith(".fasta")]

# Process each FASTA file
for fasta_file in fasta_files:
    df = read_fasta(fasta_file)
    encoded_df = df.withColumn("encoded", one_hot_udf(df["sequence"]))

    # Save image function
    def save_image(row):
        accession = row.accession.replace('.', '_')
        encoded_seq = row.encoded
        height = 4
        width = len(encoded_seq)

        # Create image array
        img_array = np.zeros((height, width), dtype=np.uint8)
        for i, nuc in enumerate(encoded_seq):
            img_array[:, i] = nuc

        img = Image.fromarray(img_array * 255, 'L')
       
        # Extract variant type from the description
        variant = row.description.split('|')[2]
        variant_dir = os.path.join("Dataset/one-hot", variant)
        os.makedirs(variant_dir, exist_ok=True)
        
        img.save(os.path.join(variant_dir, f"{accession}.png"))
        print(f"Saved {accession}.png in folder {variant_dir}")

    # Apply the save_image function to each row of the encoded DataFrame
    encoded_df.foreach(save_image)

# Stop Spark session
spark.stop()
