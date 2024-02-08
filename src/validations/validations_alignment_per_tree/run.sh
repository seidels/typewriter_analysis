#load java modules
env2lmod
module load openjdk/17.0.0_35

#enter R environment
conda activate tree_dist_env

# draw the clock rate and the insert probabilities from prior distributions
sbatch --wrap="Rscript draw_simulation_params.R"

#exit R environment
conda deactivate

#simulate trees and alignments
sbatch --job-name="tree_alignment_simulations" --array=1-100 --wrap="java -jar ~/all_beasts_2.7_2024.jar -seed \$SLURM_ARRAY_TASK_ID -version_file version.xml simulate_alignment_and_tree.xml"

#run inference
sbatch --job-name="validations" --array=1-100 --wrap="java -jar ~/all_beasts_2.7_2024.jar -seed \$SLURM_ARRAY_TASK_ID -version_file version.xml infer_given_fixed_tree_13_inserts.xml"
