#!/bin/zsh


# draw the clock rate and the insert probabilities from prior distributions
Rscript draw_simulation_params.R

# run the simulation
for seed in `seq 3 3`
do
    #java -jar 23-04-05-typewriter.jar -seed $seed -D outputDir="/Users/seidels/Projects/typewriter_analysis/results/validations/validations_alignment_per_tree/simulated_values/" simulate_alignment_and_tree.xml

    java -jar 23-04-05-typewriter.jar -overwrite -seed $seed infer_given_fixed_tree_13_inserts.xml
done
