#!/bin/bash
#SBATCH --job-name=01b
#SBATCH --output=01b.out
#SBATCH --error=01b.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=hpcf2013
#SBATCH --qos=medium
#SBATCH --mem=MaxMemPerNode

module load Keras

mpirun python ./01b_create_5p_data_test.py

