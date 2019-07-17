#!/bin/bash
#SBATCH --job-name=spring
#SBATCH --output=spring.out
#SBATCH --error=spring.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=hpcf2013
#SBATCH --qos=medium
#SBATCH --mem=MaxMemPerNode

module load Keras

mpirun python ./03a_create_spring.py

