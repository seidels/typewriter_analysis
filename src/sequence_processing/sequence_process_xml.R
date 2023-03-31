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
library(stringr)

## ---------------------------

# load the preprocessed data ---------

#edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))
#load the filtered data: 

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

#extract 100/500 cells from all targetBCs (except the one that only has 2 edits, total 12) from paper:
targetBCs = c("ATGGTAAG","ATTTATAT",
                         "ATTTGGTT", "GCAGGGTG",
                         "GTAAAGAT", "TAGATTTT",
                         "TGCGATTT", "TGGACGAC",
                         "TGGTTTTG", "TTAGATTG",
                         "TTGAGGTG",
                         "TTTCGTGA")


#here, change all_tbcs with desired vector of targetBCs.
sample_dataset_for_BEAST(100,targetBCs,edit_table_by_5,"results/analysis_cell_culture_data/simple_100cells_12tbcs/alignment.txt")
sample_dataset_for_BEAST(500,targetBCs,edit_table_by_5,"results/analysis_cell_culture_data/simple_500cells_12tbcs/alignment.txt")
sample_dataset_for_BEAST(1000,targetBCs,edit_table_by_5,"results/analysis_cell_culture_data/simple_1000cells_12tbcs/alignment.txt")



##################
#helper functions#
##################

sample_dataset_for_BEAST <- function(n_cells,targetBCs,data,name="alignment.txt") {


for(i in 1:length(targetBCs)) {
targetBC <- targetBCs[i]  
sampled <- sample_targetBCs(size=n_cells,dataset=data,targetBC)
write( paste0("<data  id=\"data_",targetBC,"\" spec=\"Alignment\" name=\"alignment\" >
            <userDataType spec=\"beast.evolution.datatype.ScarData\" nrOfStates=\"20\"/>"),name,append = TRUE)
for(j in 1:n_cells) {
  write(paste0("            <sequence spec=\"Sequence\" taxon=\"",j-1,"\"  value=\"",sampled[j],"\"/>"),name,append=TRUE)

}
write("</data>",name,append=TRUE)


}

}

# sample based on TargetBC 
# to do check: can there be several copies of the same targetBC per cell? 

sample_targetBCs <- function(size,dataset,targetBC) {
  
  #sample from sequences from this particular TargetBC
  print(targetBC)
  print(length(dataset$beast_seq[which(dataset$TargetBC == targetBC)]))
  sample <- sample(dataset$beast_seq[which(dataset$TargetBC == targetBC)],size)
  
  #return the sample
  return(sample)
  
}




#copy/paste this from the console for insert frequencies
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]),useNA = "always")) %>% arrange(desc(Freq))
bulk_insert_count <- bulk_insert_count$Freq / sum(bulk_insert_count$Freq)
cat(bulk_insert_count)