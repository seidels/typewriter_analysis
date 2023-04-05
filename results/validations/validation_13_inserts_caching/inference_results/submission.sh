#!/bin/sh
#BSUB -W 48:00
#BSUB -R "rusage[mem=4096]"
#BSUB -J "validation_caching[2-100]"
module load java
SEED=$LSB_JOBINDEX

#analysis of the complete dataset
analysis=infer_given_fixed_tree_13_inserts


java -jar 2023-03-31-typewriter.jar -overwrite -statefile ${analysis}.${SEED}.state -seed $SEED ${analysis}.xml > ${analysis}.${SEED}.out

