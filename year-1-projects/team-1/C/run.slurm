#!/bin/bash

#SBATCH --job-name=pollution                # Job name
#SBATCH --output=slurm.out                  # Output file name
#SBATCH --error=slurm.err                   # Error file name
#SBATCH --partition=batch                   # Queue (partition)
#SBATCH --nodes=1                           # Nodes
#SBATCH --ntasks-per-node=16                # MPI tasks per node requested
#SBATCH --mem=MaxMemPerNode                 # Memory per node requested
#SBATCH --constraint=hpcf2013
#SBATCH --exclusive

srun ./pollution 32 1e-9 999 24.0 1.0e-1

