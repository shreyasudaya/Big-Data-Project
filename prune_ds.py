import os
import shutil

def copy_dataset_with_cap(source_dir, target_dir, file_cap=1000):
    # Ensure target directory exists
    os.makedirs(target_dir, exist_ok=True)
    
    # Loop through each subfolder (variant) in the source directory
    for variant in os.listdir(source_dir):
        variant_path = os.path.join(source_dir, variant)
        if os.path.isdir(variant_path):
            target_variant_dir = os.path.join(target_dir, variant)
            os.makedirs(target_variant_dir, exist_ok=True)
            
            # Get all files in the variant subfolder
            csv_files = [f for f in os.listdir(variant_path) if f.endswith('.csv')]
            
            # Limit the number of files to copy
            files_to_copy = csv_files[:file_cap]
            
            # Copy the files
            for csv_file in files_to_copy:
                source_file = os.path.join(variant_path, csv_file)
                target_file = os.path.join(target_variant_dir, csv_file)
                shutil.copy2(source_file, target_file)
                print(f"Copied: {source_file} -> {target_file}")

    print("Dataset copy complete.")

# Paths
source_directory = "..\one-hot"  # Replace with your actual source path
target_directory = "Dataset\one-hot"  # Replace with your actual target path

# Copy dataset with a file cap of 1000 per variant folder
copy_dataset_with_cap(source_directory, target_directory, file_cap=1000)
