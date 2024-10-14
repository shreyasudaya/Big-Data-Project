import os
import shutil
import random

def split_dataset(source_dir, target_dir, train_ratio=0.75):
    # Create the target directory if it doesn't exist
    os.makedirs(target_dir, exist_ok=True)

    # Loop through each class folder in the source directory
    for class_folder in os.listdir(source_dir):
        class_folder_path = os.path.join(source_dir, class_folder)

        # Check if it's a directory
        if os.path.isdir(class_folder_path):
            # Create the same class folder in the target directory
            train_class_folder = os.path.join(target_dir, 'train', class_folder)
            val_class_folder = os.path.join(target_dir, 'val', class_folder)
            os.makedirs(train_class_folder, exist_ok=True)
            os.makedirs(val_class_folder, exist_ok=True)

            # List all files in the class folder
            files = os.listdir(class_folder_path)
            random.shuffle(files)  # Shuffle the files for randomness

            # Calculate the split index
            split_index = int(len(files) * train_ratio)

            # Copy files to the training folder
            for file in files[:split_index]:
                shutil.copy(os.path.join(class_folder_path, file), train_class_folder)

            # Copy files to the validation folder
            for file in files[split_index:]:
                shutil.copy(os.path.join(class_folder_path, file), val_class_folder)

if __name__ == "__main__":
    source_directory = "/home/nitk217cs010/BRICS_Viral-GNN/Dataset/Fcgr/"  # Source directory with variant classes
    target_directory = "/home/nitk217cs010/BRICS_Viral-GNN/Dataset/image-net/fcgr/"  # Target directory to create train/val splits

    split_dataset(source_directory, target_directory)
    print(f"Dataset split into training and validation sets in '{target_directory}'.")

