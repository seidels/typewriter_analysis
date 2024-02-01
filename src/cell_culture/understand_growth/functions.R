# Under a birth-death process with birth rate lambda and death rate mu, the probability
# that a lineage leaves n descendants after time t is p_n_t (Kendall 1949). The following
# functions calculate this probability


# probability to observe 0 descendants after time t.

prob_zero_descendants = function(lambda, mu, t){
  if (lambda == mu) {
    stop("lambda should not be equal to mu to avoid division by zero.")
  }

  rate_diff = mu - lambda
  numerator = mu * (1 - exp(rate_diff * t))
  denominator = lambda - mu * exp(rate_diff * t)

  return(numerator / denominator)
}

prob_one_descendant = function(lambda, mu, t){

  if (lambda == mu) {
    stop("lambda should not be equal to mu to avoid division by zero.")
  }

  rate_diff = mu - lambda

  numerator = (-rate_diff)^2 * exp(rate_diff * t)
  denominator = (lambda - mu * exp(rate_diff*t) )^2

  return(numerator / denominator)
}

prob_n_descendants = function(lambda, mu, t, n){

  if (lambda == mu) {
    stop("lambda should not be equal to mu to avoid division by zero.")
  }
  if (mu == 0) {
    stop("mu should not be 0 to avoid division by zero.")
  }


  probability = (lambda / mu * prob_zero_descendants(lambda = lambda, mu = mu, t = t))^(n-1) *
    prob_one_descendant(lambda = lambda, mu =  mu, t = t)

  return(probability)
}

## Analytical expressions based on the birth death sampling process,
## based on the analytical expressions in the doc in this folder.

prob_zero_descendants = function(lambda, mu, rho, t){

  numerator = rho * (lambda - mu)
  denominator = rho * lambda + (lambda * (1-rho) - mu) * exp(-(lambda-mu)*t)

  prob_0 = 1 - numerator/denominator

  return(prob_0)
}

prob_one_descendant = function(lambda, mu, rho, t){

  numerator = rho * (lambda - mu)^2 * exp(-(lambda-mu)*t)
  denominator = rho * lambda + (lambda * (1-rho) - mu) * exp(-(lambda-mu)*t)

  prob_1 = numerator / denominator^2

  return(prob_1)

}

q = function(lambda, mu, rho, t){

  return(prob_one_descendant(lambda, mu, rho, t) /
           (1 - prob_zero_descendants(lambda, mu, rho,t)))
}


expected_value_n = function(lambda, mu, rho, t){

  return( (1 - prob_zero_descendants(lambda, mu, rho, t)) /
            q(lambda, mu, rho, t))
}

variance_n = function(lambda, mu, rho, t){

  q_t = q(lambda, mu, rho, t)
  var_n = (1 - prob_zero_descendants(lambda, mu, rho, t)) *
    (1 - q_t) / q_t^2

  return(var_n)
}

sd_n = function(lambda, mu, rho, t){

  var_n = variance_n(lambda, mu, rho, t)

  return(sqrt(var_n))
}


## Scripts useful for manipulating trajectory data

require(tidyverse)

parseTrajectory <- function(trajStr) {
  strValues <- str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE)
  values <- apply(strValues[,-2], 2, as.numeric)
  time <- values[,1]
  src <- values[,2]
  dest <- values[,3]
  mult <- values[,4]
  N <- values[,-(1:4)]
  event<- strValues[,2]

  res <- list(time = time,
              N = N,
              event = event,
              src = src,
              dest = dest,
              mult = mult)

  return(res)
}

loadTrajectories <- function(filename, burninFrac=0.1, subsample=NA) {
  states <- NULL
  events <- NULL

  message("Loading ", filename,"...", appendLF = FALSE)
  df_in <- read_tsv(filename, col_types="ic")

  if (burninFrac>0) {
    n <- dim(df_in)[1]
    df_in <- df_in[-(1:ceiling(burninFrac*n)),]
  }

  if (!is.na(subsample)) {
    indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
    df_in <- df_in[indices,]
  }

  for (row in 1:(dim(df_in)[1])) {
    trajStr <- df_in[row,2]
    trajStates <- parseTrajectory(trajStr)
    Ndim <- dim(trajStates$N)

    if (length(Ndim)==0) {
      ntypes <- 1
      states <- bind_rows(states,
                          tibble(traj=row,
                                 type=0,
                                 time=trajStates$time,
                                 N=trajStates$N))
    } else {
      ntypes <- dim(trajStates$N)[2]
      for (s in 1:ntypes) {

        states <- bind_rows(states,
                            tibble(traj=row,
                                   type=s-1,
                                   time=trajStates$time,
                                   N=trajStates$N[,s]))
      }
    }

    events <- bind_rows(events,
                        tibble(traj=row,
                               time=trajStates$time,
                               event=trajStates$event,
                               src=trajStates$src,
                               dest=trajStates$dest,
                               mult=trajStates$mult))
  }

  states <- states %>% group_by(traj) %>% mutate(age=max(time)-time)
  events <- events %>% group_by(traj) %>% mutate(age=max(time)-time)

  message("done.")

  return(list(states=states, events=events))
}

gridTrajectories <- function(trajStates, times) {
  return(trajStates %>%
           group_by(traj, type) %>%
           summarize(N=approx(time, N, times, method="constant", f=1, yleft=0)$y,
                     age=ages))

}

gridTrajectoriesByAge <- function(trajStates, ages) {
  return(trajStates %>%
           group_by(traj, type) %>%
           summarize(N=approx(age, N, ages, method="constant", f=0, yright=0)$y,
                     age=ages))
}

