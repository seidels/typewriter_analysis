## ---------------------------
##
## Purpose of script: Implement equations 4, 8, 10 from https://arxiv.org/pdf/1401.0718.pdf
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


prob_0_events_at_time = function(mu1, t){
  prob = exp(-mu1 * t)  
  return(prob)
}

prob_1_event_at_time = function(mu1, mu2, t){
  
  prob = mu1 * (exp(-mu1 * t) / (mu2 - mu1) + exp(-mu2 * t)/ (mu1 - mu2))

  return(prob)
}

prob_2_events_at_time = function(mu1, mu2, mu3, t){
  
  prob = mu1 * mu2 * (
    
    exp(-mu1 * t)/(mu3 - mu1) / (mu2 -mu1) + 
    exp(-mu2 * t)/(mu3 -mu2) / (mu1 - mu2) +
    exp(-mu3 *t)/(mu2 - mu3)/ (mu1 - mu3)
  )
  
  return(prob)
}

