#!/bin/bash
#BSUB -W 04:00
#BSUB -R "rusage[mem=51200]"
#BSUB -J "clade_sets"
SEED=$LSB_JOBINDEX

#set environment
env2lmod
module load openjdk/17.0.0_35


~/beast/bin/applauncher CladeSetComparator -tree1 typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.tree.thinned1000000.1.trees -tree2 typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.tree.thinned1000000.2.trees -png trees1vs2.png  -verbose true


