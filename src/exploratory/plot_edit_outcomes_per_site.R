## ---------------------------
##
## Script name: plot_edit_outcomes_per_site
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

library(ggplot2)
library(reshape2)


# load the preprocessed data ---------

edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))

## ---------------------------
## explore edit outcomes per site

edits_melted = melt(edit_table_by_5, id.vars = c("Cell", "TargetBC", "nUMI"))

# plot only first 3 sites
g_3_sites =
  ggplot(data = subset(x = edits_melted, edits_melted$variable %in% c("Site1", "Site2", "Site3")), aes(x=value)) +
  facet_grid(variable ~.)+
  geom_bar()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = g_3_sites, filename = paste0(plot_path, "edit_outcomes_sites_13.jpg"))

# plot all sites for complete overview (albeit increasing loss with #sites makes it difficult
# to see edit proportions)
g_5_sites =ggplot(data = edits_melted, aes(x=value)) +
  facet_grid( variable ~.)+
  geom_bar()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = g_5_sites, filename = paste0(plot_path, "edit_outcomes_sites_15.jpg"))
