#!/bin/bash
#BSUB -W 04:00
#BSUB -R "rusage[mem=51200]"
#BSUB -J "post_processing"
SEED=$LSB_JOBINDEX

#set environment
env2lmod
module load openjdk/17.0.0_35

# run analyses
#~/beast/bin/logcombiner -log *.trees -b 10 -o typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.combined1-3.trees

#does not run, ie. gets heap space error even with 50g memory; instead use thinned trees
#~/beast/bin/treeannotator -burnin 50 -heights mean typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.combined1-3.trees typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.combined1-3.MAP.tree

~/beast/bin/logcombiner -resample 100000 -log typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.tree.?.trees -b 10 -o typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.combined1-3.thinned100000.trees

~/beast/bin/treeannotator -burnin 0 -heights mean typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.combined1-3.thinned100000.trees typewriter_partiallyPooledClock_13Sites_500Cells_DataSet1_bactrianOperator.combined1-3.thinned100000.MAP.tree



