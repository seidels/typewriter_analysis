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

setwd("/Users/seidels/Projects//typewriter_analysis/")
library(dplyr)
library(ggplot2)
library(tidyr)
library(boot)
library(cowplot)
library(stringr)

## ---------------------------

# load the preprocessed data ---------

edit_table_by_5 = readRDS("data/edit_table_filtered.RDS")


### process the dataset to a txt format for beast input: trinucleotides are mapped to an integer
### this dataset is subsampled and saved as txt file as a BEAST Alignment with ScarData


# 1st remove all trinucleotides and convert NAs with NA

for(i in 3:7) {
  edit_table_by_5[,i] <- substring(edit_table_by_5[,i], 1,3)
  edit_table_by_5[which(is.na(edit_table_by_5[,i])),i] <- "NA"

}

# 2nd, map these to an integer:

# get this from the frequencies:
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]))) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$codes <- as.character(seq(0,19))
bulk_insert_count <- data.frame(edit=bulk_insert_count$Var1, code=as.character(seq(0,19)))

pick_code <- function(edit,code_map) {
  code <- code_map[which(code_map$edit==edit),2]
  return(code)

}

# convert these trinucleotides to integer format:
for(i in 3:7) {
  edit_table_by_5[,i] <- unlist(lapply(edit_table_by_5[,i], function(x) {pick_code(x,bulk_insert_count)}))

}

# concatenate them, adding commas:
edit_table_by_5$beast_seq <- apply(edit_table_by_5[,3:7],1,function(x) {str_flatten(x,collapse = ",")})

sample_dataset_for_BEAST(100,4,edit_table_by_5, "/Users/seidels/Projects/typewriter_analysis/results/analysis_cell_culture_data/simple_100cells_4tbcs/alignment.txt")




##################
#helper functions#
##################

sample_dataset_for_BEAST <- function(n_cells,n_targetBCs,data,name="typewriter_data.txt") {

sink(name)
for(i in 1:n_targetBCs) {
sampled <- sample_targetBCs(size=n_cells,dataset=data,targetBC_index=i)
cat(paste0("<data  id=\"data_",i,"\" spec=\"Alignment\" name=\"alignment\" >
            <userDataType spec=\"beast.evolution.datatype.ScarData\" nrOfStates=\"14\"/>"))
cat("\n")
for(j in 1:n_cells) {
  cat(paste0("            <sequence spec=\"Sequence\" taxon=\"",j-1,"\"  value=\"",sampled[j],"\"/>"))
  cat("\n")
}
cat("</data>")
cat("\n")
cat("\n")

}
sink()
}

# sample based on TargetBC
# to do check: can there be several copies of the same targetBC per cell?

sample_targetBCs <- function(size,dataset,targetBC_index =1) {
  #select all sequences with a specific target BC
  targetBC <- dataset$TargetBC[targetBC_index]

  #sample from sequences from this particular TargetBC
  sample <- sample(dataset$beast_seq[which(dataset$TargetBC == targetBC)],size)

  #return the sample
  return(sample)

}





#close(fileConn)

#copy/paste this from the console for insert frequencies
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]),useNA = "always")) %>% arrange(desc(Freq))
bulk_insert_count <- bulk_insert_count$Freq / sum(bulk_insert_count$Freq)
cat(bulk_insert_count)
