#!/bin/bash

#SBATCH --job-name=ML                   # Job name
#SBATCH --output=out_err/slurm.out      # Output file name
#SBATCH --error=out_err/slurm.err       # Error file name
#SBATCH --partition=high_mem            # Partition
#SBATCH --qos=short+                    # Queue
#SBATCH --time=01:00:00                 # Time limit 
#SBATCH --nodes=8                       # Number of nodes
#SBATCH --ntasks-per-node=8             # MPI processes per node

mpiexec python par_cnn.py
