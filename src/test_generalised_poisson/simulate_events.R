## ---------------------------
##
## Purpose of script: Check the correctness of equation 10 of the generalised
## poisson process from https://arxiv.org/pdf/1401.0718.pdf
##
## Author: Sophie Seidel
##
## Date Created: 2024-03-06
##
## Copyright (c) Sophie Seidel, 2024
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac

setwd("~/Projects/typewriter_analysis/src/test_generalised_poisson/")      # Sophie's working directory (mac)

## ---------------------------

event_rates = c(1, 2, 3, 0.5, 1.5, 0)

n_trials = 100000
n_end = 5
t_end = 10

t = 0
n = 1

arrival_time_dat = data.frame(t1 = NA, t2= NA, t3 = NA, t4 = NA, t5 = NA)

# draw time until the next event from exponential distribution
draw_arrival_times = function(event_rates, t_end, n_end){
  
  arrival_times = rep(NA, 5)
  
  t=0
  n=1
  
  while (t < t_end & n <= n_end) {
    
    t_next_event = rexp(n = 1, event_rates[n])
    
    t = t + t_next_event
    arrival_times[n] = t
    
    #print(n)
    #print(t)
    
    n = n + 1
    
  }
 
  return(arrival_times) 
}

wt = draw_arrival_times(event_rates, t_end = 10, n_end = 5)
wt

for(trial in 1:n_trials){
  
  arrival_times = draw_arrival_times(event_rates, t_end = 10, n_end = 5)
  arrival_time_dat[trial, ] = arrival_times
}

## ---------------------------
# Arrivals per time unit

n_events_at_time_dat = data.frame(matrix(nrow = n_trials, ncol = 10))

get_n_events_at_time = function(arrival_times){
  
  n_events_at_time = rep(NA, 10)

    for (time in 1:10){
      n_events = max(which(arrival_times < time))
      n_events_at_time[time] = n_events
    }
  
  return(n_events_at_time)
}

for(trial in 1:n_trials){
  
  arrival_times = arrival_time_dat[trial, ]
  n_events_at_time = get_n_events_at_time(arrival_times)
  n_events_at_time[is.infinite(n_events_at_time)] = 0 # Correct to 0
  
  n_events_at_time_dat[trial, ] = n_events_at_time
}
## ---------------------------

table(n_events_at_time_dat[,1])/n_trials
prob_0_events_at_time(mu1 = 1, t = 1)
prob_1_event_at_time(mu1 = 1, mu2 = 2, 1)
prob_2_events_at_time(mu1 = 1, mu2 = 2, mu3 = 3, t = 1)

table(n_events_at_time_dat[,3])/n_trials
prob_0_events_at_time(mu1 = 1, t = 3)
prob_1_event_at_time(mu1 = 1, mu2 = 2, 3)
prob_2_events_at_time(mu1 = 1, mu2 = 2, mu3 = 3, t = 3)



