#!/bin/bash
#SBATCH --job-name=5p_2d32v
#SBATCH --output=5p_2d32v.out
#SBATCH --error=5p_2d32v.err
#SBATCH --partition=gpu
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --mem=MaxMemPerNode
#SBATCH --gres=gpu:2

# I'm not going to deal with conflicts. Destroy them.
module purge
module load shared default-environment

# Naiively setup andaconda as if it was 2nd python
export CONDA_PREFIX=/umbc/xfs1/cybertrn/cybertraining2019/team5/Research/miniconda3
export PYTHONPATH=$CONDA_PREFIX/bin:$PYTHONPATH
export PATH=$CONDA_PREFIX/bin:$PATH

# The CUDA toolkit location
TOOLKITDIR=/umbc/xfs1/cybertrn/common/nvidia-toolkit/
# The end-user drivers for the VXXX's GPUs
DRIVERDIR=$TOOLKITDIR
# Establish the environment
export PATH="$PATH:$TOOLKITDIR/bin"
export CUDA_HOME="$TOOLKITDIR" && \
export CUDA_PATH="$TOOLKITDIR" && \
export CUDA_ROOT="$TOOLKITDIR" && \
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOOLKITDIR/lib64" && \
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOOLKITDIR/extras/CUPTI/lib64" && \
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOOLKITDIR/nvvm/lib64" && \
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$DRIVERDIR/NVIDIA-Linux-x86_64-396.44/" && \
export LIBRARY_PATH="$LIBRARY_PATH:$TOOLKITDIR/lib64" && \
export LIBRARY_PATH="$LIBRARY_PATH:$TOOLKITDIR/lib64/stubs" && \
export PATH="$PATH:$TOOLKITDIR/bin" && \
export PATH="$PATH:$TOOLKITDIR/nvvm/bin"


hostname
nvcc --version
mpirun python ./5p_2d32v.py
