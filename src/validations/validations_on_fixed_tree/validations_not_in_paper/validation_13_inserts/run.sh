#!/bin/zsh

for seed in `seq 1 100`	    
do
    java -jar $HOME/frameworks/new_beasts/out/artifacts/organoid_jar/organoid.jar -seed $seed simulate_alignment_given_tree.xml
done
