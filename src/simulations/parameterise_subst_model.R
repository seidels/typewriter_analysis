## ---------------------------
##
## Script name: get_edit_proportions
##
## Purpose of script: Plot the number of edits introduced at different sites
## based on data in https://doi.org/10.1038/s41586-022-04922-8
##
## Author: Sophie Seidel
##
## Date Created: 2022-08-22
##
## Copyright (c) Sophie Seidel, 2022
## Email: sophie.seidel@posteo.de
##
## ---------------------------
##
## Notes: data obtained from https://github.com/shendurelab/DNATickerTape/
##
##
## ---------------------------

## set working directory for Mac
setwd("~/Projects/typewriter_analysis/")
library(ggplot2)
library(reshape2)


# load the preprocessed data ---------
edit_table_filtered = readRDS(file = "data/edit_table_filtered.RDS")


## --------
edit_outcomes_site_1 = as.data.frame(table(as.character(edit_table_filtered$Site1), useNA = "ifany"))

nr_of_edited = sum(edit_outcomes_site_1[which(!(is.na(edit_outcomes_site_1$Var1))), "Freq"])
nr_of_unedited = edit_outcomes_site_1[which(is.na(edit_outcomes_site_1$Var1)), "Freq"]
nr_of_all_sites = nr_of_edited + nr_of_unedited

fraction_edited = nr_of_edited / nr_of_all_sites

edit_duration = 25

get_rate_from_fraction_edited(fraction_edited, edit_duration)

mean_number_edits_for_25d = get_rate_from_fraction_edited(fraction_edited, edit_duration) * 25


edit_table_filtered$nrOfEdits = get_nr_of_edits(edit_table_filtered)

ggplot(edit_table_filtered, aes(x=nrOfEdits)) +
  geom_histogram()
## NOTE: This histogram looks like it is the overlay of several distributions

## Idea: separate out TargetBCs that show limited editing
target_barcodes = unique(edit_table_filtered$TargetBC)
edits_per_target_barcode = data.frame(target_barcode = target_barcodes, max_edits = 0)

for (target_barcode in target_barcodes){
  print(target_barcode)
  target_barcode_data = edit_table_filtered[which(edit_table_filtered$TargetBC == target_barcode), ]
  edits_per_target_barcode[which(edits_per_target_barcode$target_barcode == target_barcode), "max_edits"] = max(target_barcode_data$nrOfEdits)
}

## This is only valid for target barcode:
edits_per_target_barcode[which(edits_per_target_barcode$max_edits < 5), ]
##... which has at max 3 introductions

## Hence split the histogram with the number of edits per target bc
g = ggplot(edit_table_filtered, aes(x=nrOfEdits)) +
  facet_wrap(~ TargetBC)+
  geom_histogram()+
  theme_minimal()
g
## TODO Control for the number of cells by duplication
edit_table_filtered_reduced = unique(edit_table_filtered[, 2:7])
edit_table_filtered_reduced$nrOfEdits = get_nr_of_edits(edit_table_filtered_reduced)

g = ggplot(edit_table_filtered_reduced, aes(x=nrOfEdits)) +
  facet_wrap(~ TargetBC)+
  geom_histogram()+
  theme_minimal()

get_rate_from_fraction_edited = function(fraction_edited, edit_duration){

  edit_introduction_rate = - log(1 - fraction_edited) / edit_duration

  return(edit_introduction_rate)
}

get_nr_of_edits = function(edit_table){
  ## compare to data
  nr_of_edits = rep(x = 0, nrow(edit_table))

  for (i in 1:nrow(edit_table)){

    if(is.na(edit_table[i, "Site1"])){

      nr_of_edits[i] = 0
      next

    }else if(is.na(edit_table[i, "Site2"])){

      nr_of_edits[i] = 1
      next

    }else if (is.na(edit_table[i, "Site3"])){

      nr_of_edits[i] = 2
      next

    }else if (is.na(edit_table[i, "Site4"])){

      nr_of_edits[i] = 3
      next

    }else if (is.na(edit_table[i, "Site5"])){

      nr_of_edits[i] = 4
      next

    }else{
      nr_of_edits[i] = 5
    }
  }
  return(nr_of_edits)
}
