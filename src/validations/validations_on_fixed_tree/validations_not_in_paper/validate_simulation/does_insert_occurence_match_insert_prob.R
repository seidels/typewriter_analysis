## ---------------------------
##
## Script name: validate simulation
##
## Purpose of script:
##
## Author: Sophie Seidel
##
## Date Created: 2023-01-11
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

setwd("~/Projects/")      # Sophie's working directory (mac)

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R")

## ---------------------------
library(stringr)

get_data_from_child_sequences = function(child_1, child_2){

  child_1_data = sapply(X = c("1", "2"), function(x){
    str_count(string = child_1, pattern = x)
  } )
  child_2_data = sapply(X = c("1", "2"), function(x){
    str_count(string = child_2, pattern = x)
  } )

  data = data.frame(child_1 = child_1_data, child_2=child_2_data )

  return(data)
}

### simulation 1
#### clock rate = 2; insert probs = c(0.5 0.5)
child_1 = "1,1,2,1,1,1,1,2,2,2,1,2,1,2,1,2,1,2,1,1,1,1,2,2,2,1,1,2,2,1,1,1,2,2,1,2,1,1,1,2,2,2,1,2,1,1,2,2,2,1,1,2,0,0,0"
child_2 = "2,1,1,1,1,2,1,2,1,2,1,2,1,2,1,2,2,1,2,2,1,1,2,2,1,2,2,1,2,1,2,1,2,2,1,1,1,1,2,2,1,2,2,1,1,1,1,1,1,2,1,1,2,2,2"

data = get_data_from_child_sequences(child_1, child_2)

insert_occs = rowSums(data)
insert_occs / sum(insert_occs)

### simulation 2
#### clock rate=4; insert probs = c(0.5, 0.5)
child_1 = "1,1,2,1,1,1,1,2,2,2,1,2,1,2,1,2,1,2,1,1,1,1,2,2,2,1,1,2,2,1,1,1,2,2,1,2,1,1,1,2,2,2,1,2,1,1,2,2,2,1,1,2,1,2,2,1,1,1,1,2,1,2,1,2,1,2,1,2,1,2,2,1,2,2,1,1,2,2,1,2,2,1,2,1,2,1,2,2,1,1,1,1,2,2,1,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
child_2 = "1,1,2,2,2,2,1,1,2,2,1,2,2,2,1,1,2,1,1,1,1,2,2,2,1,2,2,2,1,2,1,1,1,2,2,1,1,2,1,2,1,1,2,1,2,1,2,2,2,2,1,2,2,1,2,2,2,2,2,1,1,1,2,1,1,2,2,1,1,1,1,1,2,1,2,2,2,1,2,1,2,1,2,1,2,1,1,1,1,2,2,2,2,2,1,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"

data = get_data_from_child_sequences(child_1, child_2)

insert_occs = rowSums(data)
insert_occs / sum(insert_occs)

sum(data)
4*24*2

### simulation 2
#### clock rate=4; insert probs = c(0.1, 0.9)

child_1 = "2,2,2,2,2,1,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,1,2,2,2,2,2,2,2,1,2,1,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
child_2 = "2,1,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,2,2,2,2,2,1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,1,2,2,2,2,2,1,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,1,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"

data = get_data_from_child_sequences(child_1, child_2)

insert_occs = rowSums(data)
insert_occs / sum(insert_occs)

sum(data)
4 * 24 * 2
