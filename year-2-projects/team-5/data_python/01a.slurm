#!/bin/bash
#SBATCH --job-name=01a
#SBATCH --output=01a.out
#SBATCH --error=01a.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=hpcf2013
#SBATCH --qos=medium
#SBATCH --mem=MaxMemPerNode

module load Keras

mpirun python ./01a_create_5p_data_train.py

