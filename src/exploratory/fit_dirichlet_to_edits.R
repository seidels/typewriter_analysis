## ---------------------------
##
## Script name: fit_dirichlet_to_edits
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
library(compositions)

## set working directory for Mac
setwd("~/Projects/typewriter_analysis/")
library(ggplot2)
library(reshape2)


# load the preprocessed data ---------
edit_table_filtered = readRDS(file = "data/edit_table_filtered_reduced.RDS")
frequent_insertions = read.csv(file = "src/exploratory/13_frequent_TargetBCs.csv")


edit_outcomes_site_1 = as.data.frame(table(as.character(edit_table_filtered$Site1), useNA = "ifany"))

target_bcs = unique(edit_table_filtered$TargetBC)
insert_bcs = unique(edit_table_filtered$Site1)
insert_bcs = insert_bcs[!is.na(insert_bcs)]

n_target_bcs = length(target_bcs)
n_insert_bcs = length(insert_bcs)

compositions = matrix(nrow = n_target_bcs, ncol = n_insert_bcs, data = 0 )

for (i in 1:n_target_bcs){

  inserts_for_target_bc = edit_table_filtered[which(edit_table_filtered$TargetBC == target_bcs[i]), "Site1"]

  insert_occurences = as.data.frame(table(inserts_for_target_bc))

  #subset for inserts
  insert_composition = data.frame(insert_bcs = insert_bcs, occurences = insert_occurences[match(insert_bcs, insert_occurences$inserts_for_target_bc), "Freq"])
  insert_composition[is.na(insert_composition$occurences), "occurences"] = 0
  insert_composition$normalised = insert_composition$occurences / sum(insert_composition$occurences)

  compositions[i, ] = insert_composition$normalised
}

comp_df = as.data.frame(compositions)
colnames(comp_df) = insert_bcs
comp_df$target_bc = target_bcs

library(reshape2)
comp_melt = melt(comp_df,id.vars = "target_bc")
ggplot(comp_melt, aes(x=variable, y=value)) +
  geom_col()+
  facet_grid(target_bc~.)+
  xlab("insert bc") +
  ylab("proportion")

# fit dirichlet to composition of edits
library(Compositional)
## add super small constant s.t. 0 entries are non-zero
fits = diri.est(x = (compositions +0.0000001), type = "prec")

rDirichlet.acomp(n = 1,  alpha = fitted$alpha)
