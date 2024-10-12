#!/bin/sh
#SBATCH --job-name=obtain_img          # Job name
#SBATCH --partition=gpu                   # Partition name (GPU partition)
#SBATCH --gres=gpu:2                      # Number of GPUs required
#SBATCH --nodes=1                         # Number of compute nodes
#SBATCH --ntasks=2                        # Number of tasks
#SBATCH --cpus-per-task=4                 # Number of CPU cores per task
#SBATCH --mem=16G                         # Total memory per node
#SBATCH --time=24:00:00                   # Time limit (HH:MM:SS)
#SBATCH --output=kmeronehot.out               # Standard output log
#SBATCH --error=kmeronehot.err                # Standard error log

# Load necessary modules
module load cuda/11.7
module load anaconda3/pytorch
module load horovod_python/3.9
module load openmpi/4.1.4 
# Activate Conda environment if needed


# Upgrade pip to the latest version
#pip install --upgrade pip

# Install required Python packages
pip install --user numpy
pip install --user biopython    # Assuming 'Bio' refers to Biopython
pip install --user matplotlib
pip install --user futures   
pip install --user numba
pip install --user seaborn   
# Verify installations (optional)
echo "Installed Python packages:"
pip list

# Ensure that Python uses the GPU
# This depends on your scripts; ensure they utilize CUDA-enabled libraries like PyTorch

# Run your Python scripts
CUDA_VISIBLE_DEVICES=0 python proc.py &
CUDA_VISIBLE_DEVICES=1 python kprocheatmap.py

