#!/bin/bash
#SBATCH --job-name=fall
#SBATCH --output=fall.out
#SBATCH --error=fall.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=hpcf2013
#SBATCH --qos=medium
#SBATCH --mem=MaxMemPerNode

module load Keras

mpirun python ./03c_create_fall.py

