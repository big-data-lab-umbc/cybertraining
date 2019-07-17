#!/bin/bash
#SBATCH --job-name=02a
#SBATCH --output=02a.out
#SBATCH --error=02a.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=hpcf2013
#SBATCH --qos=medium
#SBATCH --mem=MaxMemPerNode

module load Keras

mpirun python ./02a_create_20p_data_train.py

