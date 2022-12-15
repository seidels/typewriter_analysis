/## ---------------------------
##
## Script name: sequence_process_xml
##
## Purpose of script: Pre-processing and sampling of sequences for xml inputs
##
## Author: Antoine Zwaans
##
## Date Created: 2022-12-13
##
## Copyright (c) Antoine Zwaans, 2022
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("/Users/azwaans/typewriter_analysis/")    
library(dplyr) 
library(ggplot2)
library(tidyr)
library(boot)
library(cowplot)

## ---------------------------

# load the preprocessed data ---------

edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))

# 1st remove all trinucleotides and convert NAs with NA

for(i in 3:7) {
  edit_table_by_5[,i] <- substring(edit_table_by_5[,i], 1,3)
  edit_table_by_5[which(is.na(edit_table_by_5[,i])),i] <- "NA"
  
}



# 2nd, map these to an integer: 

# get this from the frequencies: 
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]),useNA = "always")) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$codes <- as.character(seq(0,20))
bulk_insert_count <- data.frame(bulk_insert_count$Var1, code=as.character(seq(0,20)))

pick_code <- function(edit,code_map) {
  code <- code_map[which(code_map$codes==edit),2]
  return(code)
  
}

# convert these trinucleotides to integer format:
for(i in 3:7) {
  edit_table_by_5[,i] <- unlist(lapply(edit_table_by_5[,i], function(x) {pick_code(x,bulk_insert_count)}))
  
}


# concatenate them, adding commas: 
edit_table_by_5$beast_seq <- apply(edit_table_by_5[,3:7],1,function(x) {str_flatten(x,collapse = ",")})


# sample based on TargetBC 
# can there be several copies of the same targetBC per cell? 

sample_dataset_for_BEAST <- function(size,dataset,targetBC_index =1) {
  #select all sequences with a specific target BC
  targetBC <- dataset$TargetBC[targetBC_index]
  
  #sample from sequences from this particular TargetBC
  sample <- sample(dataset$beast_seq[which(dataset$TargetBC == targetBC)],size)
  
  #return the sample
  return(sample)
  
}

########

#copy/paste the following output from the console for sequences in xml
# here we sample 30 sequences from the 1st TargetBC available (default)


sampled <- sample_dataset_for_BEAST(30,edit_table_by_5)
for(i in 1:30) {
  cat(paste0("<sequence id=\"",i-1,"\" spec=\"Sequence\" taxon=\"",i-1,"\"  value=\"",sampled[i],"\"/>"))
  cat("\n")
}

#copy/paste this from the console for insert frequencies
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]),useNA = "always")) %>% arrange(desc(Freq))
bulk_insert_count <- bulk_insert_count$Freq / sum(bulk_insert_count$Freq)
cat(bulk_insert_count)