## ---------------------------
##
## Script name: get summary data from all log files
##
## Purpose of script:
##
## Author: Sophie Seidel
##
## Date Created: 2022-12-21
##
## Copyright (c) Sophie Seidel, 2022
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

log_dir = "./results/validation_2_inserts/inference_results/"

simulation_dir = "results/validation_2_inserts/simulation_parameters/"

parameters_of_interest = c("clockRate", "insertRate.1")


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
insert_rate_inference = data.frame(seed=1:100, hpd_lower=0, hpd_upper=0, median=0, true_value=0, recovered=F)


for (seed in  1:100){

  # get inference log
  log_file = paste0("infer_given_fixed_tree_2_inserts.", seed, ".log")
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
  true_clock_rate = simulation_parameters[3, 2]
  recovered = is_parameter_in_hpd(hpd_lower = clock_rate_hpd["lower"], hpd_upper = clock_rate_hpd["upper"], 
                      true_parameter = true_clock_rate)
  clock_rate_inference[seed, ] = c(seed, clock_rate_hpd, median(log_data_wo_burnin$clockRate), true_clock_rate, recovered)
  
  insert_rate_hpd = hdi(object = log_data_wo_burnin$insertRates.1, credMass = 0.95)
  true_insert_rate = simulation_parameters[1, 2]
  recovered = is_parameter_in_hpd(hpd_lower = insert_rate_hpd["lower"], hpd_upper = insert_rate_hpd["upper"], 
                                  true_parameter = true_insert_rate)
  insert_rate_inference[seed, ] = c(seed, insert_rate_hpd, median(log_data_wo_burnin$insertRates.1), true_insert_rate, recovered)
  
}

clock_rate_inference = clock_rate_inference[order(clock_rate_inference$true_value), ]
clock_rate_inference$orderedSeed = 1:100

g = ggplot(clock_rate_inference, aes(x=orderedSeed, y=median))+
  geom_point()+
  geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_point(aes(x=orderedSeed, y=true_value), col="darkgreen")+
  theme_classic()+
  xlab("Simulation seeds ordered with increasing simulation parameter") + 
  ylab("Estimated posterior intervals and medians")

g
ggsave(plot = g, "./src/validation_2_inserts/inference_clock_rate.pdf")

insert_rate_inference = insert_rate_inference[order(insert_rate_inference$true_value), ]
insert_rate_inference$orderedSeed = 1:100

g = ggplot(insert_rate_inference, aes(x=orderedSeed, y=median))+
  geom_point()+
  geom_errorbar(aes(ymin = hpd_lower, ymax=hpd_upper), alpha=0.4)+
  geom_point(aes(x=orderedSeed, y=true_value), col="darkgreen")+
  theme_classic()+
  xlab("Simulation seeds ordered with increasing simulation parameter") + 
  ylab("Estimated posterior intervals and medians")

g
ggsave(plot = g, "./src/validation_2_inserts/inference_insert_rate.pdf")

sum(insert_rate_inference$recovered)
sum(clock_rate_inference$recovered)
