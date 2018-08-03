#!/bin/bash                 
#SBATCH --job-name=cloud3d
#SBATCH -N 1
#SBATCH --partition=batch	  
#SBATCH --qos=medium        
#SBATCH --mem-per-cpu=5000

srun ./submit.sh 
