#!/bin/bash
#BSUB -W 120:00
#BSUB -R "rusage[mem=6144]"
#BSUB -n 5
#BSUB -J "cell_culture_500_cells_seed_[1-3]"
SEED=$LSB_JOBINDEX

#set environment
env2lmod
module load openjdk/17.0.0_35

# run analyses

java -jar ~/all_beasts2.7.jar -seed $SEED -threads 5 -overwrite typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.xml

