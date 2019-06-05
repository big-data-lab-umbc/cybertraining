#!/bin/bash

#SBATCH --job-name=test_python_team_2   # Job name
#SBATCH --output=slurm.out              # Output file name
#SBATCH --error=slurm.err               # Error file name
#SBATCH --partition=develop             # Partition
#SBATCH --qos=short                     # Queue
#SBATCH --time=01:00:00                 # Time limit 
#SBATCH --mem=4000
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # MPI processes per node

#time main.py data_eof_3.csv
#time ../main.py data_eof_5.csv
#time ../main.py data_eof_30.csv
python main.py Data_Selected_Modes.csv
time main.py Data_Selected_Modes.csv
python post_process.py
