## ---------------------------
##
## Script name: get summary data from all log files
##
## Purpose of script: Analyse well-calibrated simulations
##
## Author: Sophie Seidel
##
## Date Created: 2023-03-31
##
## Copyright (c) Sophie Seidel, 2023
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## set working directory for Mac

setwd("~/Projects/typewriter_analysis/")
## ---------------------------


## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tracerer)
library(HDInterval)
library(ggplot2)

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R")

## ---------------------------

log_dir = "/Volumes/stadler/People/Sophie_Antoine_shared/validations/validation_13_inserts_caching/inference_results/"

simulation_dir = "/Volumes/stadler/People/Sophie_Antoine_shared/validations/validation_13_inserts_caching/simulation_parameters/"

parameters_of_interest = c("clockRate", paste0("insertRate.", 1:13))


inc = function(x){
  eval.parent(substitute(x <- x+1))
}

is_parameter_in_hpd = function(hpd_lower, hpd_upper, true_parameter){
  if(true_parameter >= hpd_lower && true_parameter <= hpd_upper)
    return(TRUE)
  else
    return(FALSE)
}

nr_converged_chains = 0

clock_rate_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)
insert_rate_inference = data.frame(seed=rep(1:100, each = 13), hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F,
                                   insertRate= rep(1:13,100))


for (seed in  1:100){

  print(seed)
  # get inference log
  log_file = paste0("infer_given_fixed_tree_13_inserts.", seed, ".log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.1)

  # get true parameter from simulation
  simulation_parameter_file = paste0(simulation_dir, "simParams_", seed, ".csv")
  simulation_parameters = read.csv(simulation_parameter_file)
  colnames(simulation_parameters) = c("parameter", "value")

  #check that chains has converged
  esses = calc_esses(log_data_wo_burnin, sample_interval = 1000)
  esses = esses[!colnames(esses) %in% c("birthRate", "deathRate", "prior")]
  esses = esses[! names(esses) %in% c("samplingProportion", "treeHeight.t.alignment", "treeLength.t.alignment")]
  esses

  if (any(esses < 200)){
    problematic_ess = esses[which(esses < 200)]
    print(paste("Ess for log file ", log_file, " has ESS < 200: for parameters", paste(names(problematic_ess), collapse = ",")))
    next

  }else{
    inc(nr_converged_chains)
  }

  clock_rate_hpd = hdi(object = log_data_wo_burnin$clockRate, credMass = 0.95)
  true_clock_rate = simulation_parameters[14, 2]
  recovered = is_parameter_in_hpd(hpd_lower = clock_rate_hpd["lower"], hpd_upper = clock_rate_hpd["upper"],
                      true_parameter = true_clock_rate)
  clock_rate_inference[seed, ] = c(seed, clock_rate_hpd, median(log_data_wo_burnin$clockRate), true_clock_rate, recovered)

  for (insert_rate_nr in 1:13){

    insert_rate = paste0("insertRates.", insert_rate_nr)
    insert_rate_hpd = hdi(object = log_data_wo_burnin[, insert_rate], credMass = 0.95)
    median_insert_rate = median(log_data_wo_burnin[, insert_rate])
    true_insert_rate = simulation_parameters[insert_rate_nr, "value"]
    recovered = is_parameter_in_hpd(hpd_lower = insert_rate_hpd["lower"], hpd_upper = insert_rate_hpd["upper"],
                                    true_parameter = true_insert_rate)


    row_index = which(insert_rate_inference$seed == seed & insert_rate_inference$insertRate == insert_rate_nr)
    insert_rate_inference[row_index, ] = c(seed, as.numeric(insert_rate_hpd), median_insert_rate, true_insert_rate,
                                           recovered, insert_rate_nr)
  }
}

clock_rate_inference = clock_rate_inference[order(clock_rate_inference$true_value), ]
clock_rate_inference$orderedSeed = 1:100

saveRDS(object = clock_rate_inference, file = "/Volumes/stadler/People/Sophie_Antoine_shared/validations/validation_13_inserts_caching/inference_results/summary_data_clock.RDS")
saveRDS(object = insert_rate_inference, file = "/Volumes/stadler/People/Sophie_Antoine_shared/validations/validation_13_inserts_caching/inference_results/summary_data_insert.RDS")

g = ggplot(clock_rate_inference, aes(x=orderedSeed, y=median))+
  geom_point()+
  geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_point(aes(x=orderedSeed, y=true_value), col="darkgreen")+
  theme_classic()+
  xlab("Simulation seeds ordered with increasing simulation parameter") +
  ylab("Estimated posterior intervals and medians")

g
ggsave(plot = g, "./src/validation_13_inserts_caching/inference_clock_rate.pdf")


g = ggplot(insert_rate_inference, aes(x=seed, y=median))+
  geom_point()+
  facet_grid(insertRate~.)+
  geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_point(aes(x=seed, y=true_value), col="darkgreen")+
  theme_classic() +
  xlab("Simulation seeds") +
  ylab("Estimated posterior intervals and medians")
g

ggsave(plot = g, "./src/validation_13_inserts_caching/inference_insert_rate.pdf", width = 20, height = 30, units = "cm")

#coverages
for (insert in 1:13){
  print(sum(insert_rate_inference[insert_rate_inference$insertRate == insert, "recovered"]))
}
sum(clock_rate_inference$recovered)/nr_converged_chains

# correlation
cor.test(clock_rate_inference$median, clock_rate_inference$true_value)
cor.test(insert_rate_inference$median, insert_rate_inference$true_value)


g = ggplot(insert_rate_inference, aes(x=true_value, y=median)) +
  geom_point()+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median and posterior interval")+
  theme_classic()+
  xlim(0, 0.5) + ylim(0, 0.5)
g
ggsave(plot = g, "./src/validation_13_inserts_caching/inference_insert_rate_v2.pdf",
       width = 12, height = 12, units = "cm")

g = ggplot(clock_rate_inference, aes(x=true_value, y=median)) +
  geom_point()+geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_abline(slope = 1, col="darkgreen")+
  xlab("True value")+
  ylab("Estimated median and posterior interval")+
  theme_classic()+
  xlim(0, 0.5) + ylim(0, 0.5)
g
ggsave(plot = g, "./src/validation_13_inserts_caching/inference_clock_rate_v2.pdf",
       width = 12, height = 12, units = "cm")

