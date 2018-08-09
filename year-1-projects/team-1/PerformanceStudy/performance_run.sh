#!/bin/bash
function run_script
{
    STUDY_NAME=$(printf 'N%s' ${pollution})
    DIR_NAME=$( printf "%s/n%04dppn%04d" ${STUDY_NAME} $nodes $ppn )
    cd $DIR_NAME
    sbatch run.slurm
    cd $CWD
}

CWD=$(pwd)
for nodes in 1 2 4 8 16
do
 for ppn in 1 2 4 8 16
 do
   for N in 31 63 127 255 2047 4095 8191
   do
     pollution=$( printf "%d" $N )
     run_script
   done
 done
done
