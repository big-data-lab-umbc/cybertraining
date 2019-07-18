#!/bin/bash
#SBATCH --job-name=02b
#SBATCH --output=02b.out
#SBATCH --error=02b.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=hpcf2013
#SBATCH --qos=medium
#SBATCH --mem=MaxMemPerNode

module load Keras

mpirun python ./02b_create_20p_data_test.py

