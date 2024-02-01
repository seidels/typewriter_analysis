 
 sbatch --mem-per-cpu=40G --time=4:00:00 --wrap="java -Xms20g -Xmx40g -jar $HOME/beasts2.7.jar -overwrite -seed 1 -D nParticles=100  traj_per_state.xml"
 
sbatch --mem-per-cpu=70G --time=2:00:00 --wrap="Rscript bdmm_prime_trajectories.R"
 
 sbatch --mem-per-cpu=40G --time=4:00:00 --wrap="java -Xms20g -Xmx40g -jar $HOME/beasts2.7.jar -overwrite -seed 1 -D nParticles=100  traj_per_state_bdsky.xml"
