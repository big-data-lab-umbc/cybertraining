#!/bin/bash

#SBATCH --job-name=4kdata               # Job name
#SBATCH --output=slurm.out              # Output file name
#SBATCH --error=slurm.err               # Error file name
#SBATCH --partition=high_mem            # Partition
#SBATCH --qos=short+                    # Queue
#SBATCH --time=01:00:00                 # Time limit 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # MPI processes per node

python test_RT_wrap.py
