#!/bin/bash

#SBATCH --job-name=MOD-CAL-6m
#SBATCH --output=slurm-MOD-CAL-6m.out
#SBATCH --partition=batch
#SBATCH --qos=medium
#SBATCH --mem=16000


echo ${job} $HOSTNAME
python spot_MODIS_3Deffects.py 6 30
