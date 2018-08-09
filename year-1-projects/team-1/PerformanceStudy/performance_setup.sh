#!/bin/bash
function generate_script
{
    STUDY_NAME=$(printf 'N%s' ${pollution})
    DIR_NAME=$( printf "%s/n%04dppn%04d" ${STUDY_NAME} $nodes $ppn )
    echo "Creating job $DIR_NAME"
    mkdir -p $DIR_NAME
    
    cat << _EOF_ > ${DIR_NAME}/run.slurm
#!/bin/bash

#SBATCH --job-name=PollutionDiffusion       # Job name
#SBATCH --output=slurm.out                  # Output file name
#SBATCH --error=slurm.err                   # Error file name
#SBATCH --partition=batch                   # Queue (partition)
#SBATCH --nodes=${nodes}                    # Nodes
#SBATCH --ntasks-per-node=${ppn}            # MPI tasks per node requested
#SBATCH --mem=MaxMemPerNode                 # Memory per node requested
#SBATCH --constraint=hpcf2013
#SBATCH --exclusive

srun ../../pollution ${pollution} 1e-9 999 24.0 1.0e-1

_EOF_
}

for nodes in 1 2 4 8 16
do
    for ppn in 1 2 4 8 16
    do
        #for pow in 5 6 7 8 9 10 11   
        for N in 31 63 127 255 2047 4095 8191 
        do
        pollution=$(printf "%d" $N )
	generate_script
        done
    done
done
