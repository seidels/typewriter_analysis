#!/bin/zsh

for seed in `seq 1 100`	    
do
    java -jar 2023-03-31-typewriter.jar -seed $seed simulate_alignment_given_tree.xml
done
