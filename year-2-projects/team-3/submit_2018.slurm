#!/bin/bash

#SBATCH --job-name=Study_2018
#SBATCH --output=slurm_2018.out
#SBATCH --error=slurm_2018.err
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --nodes=1

srun python3 run_2018.py
