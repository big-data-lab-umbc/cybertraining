#!/bin/bash
#SBATCH --job-name=cloudClassifyKMeans
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=batch
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=MaxMemPerNode
#SBATCH --constraint=hpcf2013
#SBATCH --exclusive

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=16

srun python run_k-means_py3.py
