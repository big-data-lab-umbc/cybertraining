#!/bin/bash 
#SBATCH --job-name=dust_detection_parallel
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=batch
#SBATCH --nodes=4                       
#SBATCH --ntasks-per-node=16
#SBATCH --mem=MaxMemPerNode
#SBATCH --exclusive

python Deeplearning_dust_detection_feature_selection_parallel.py

