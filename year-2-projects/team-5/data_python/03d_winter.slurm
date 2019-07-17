#!/bin/bash
#SBATCH --job-name=winter
#SBATCH --output=winter.out
#SBATCH --error=winter.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=hpcf2013
#SBATCH --qos=medium
#SBATCH --mem=MaxMemPerNode

module load Keras

mpirun python ./03d_create_winter.py

