#!/bin/zsh

for seed in `seq 1 100`
do
    java -jar 23-04-05-typewriter.jar -seed $seed simulate_alignment_given_tree.xml
    java -jar 23-04-05-typewriter.jar -seed $seed simulate_alignment_length2_given_tree.xml
done
