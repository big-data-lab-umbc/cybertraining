#!/bin/bash

#SBATCH --job-name=Study_2013
#SBATCH --output=slurm_2013.out
#SBATCH --error=slurm_2013.err
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --nodes=1

srun python3 run_2013.py
