## ---------------------------
##
## Script name: plot_edit_outcomes_per_site
##
## Purpose of script: Plot the number of edits introduced at different sites
## based on data in https://doi.org/10.1038/s41586-022-04922-8
##
## Author: Sophie Seidel & Antoine Zwaans
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


## ---------------------------

## Joint counts of edits at site 1 and 2
g_count_2_and_1 = ggplot(data = edit_table_by_5)+
  facet_grid(Site1 ~.) +
  geom_bar(aes(x=Site2))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
g_count_2_and_1
ggsave(plot = g_count_2_and_1, filename = paste0(plot_path, "count_edit_2_and_1.jpg"))

## Conditional probabilities

### site 2 given 1
p_2_given_1 = conditional_a_given_b(edit_table = edit_table_by_5, site_number_a = 2, site_number_b = 1)
g_cond_21 = plot_conditional_a_given_b(conditional = p_2_given_1, site_a = "Site2", site_b = "Site1")
g_cond_21
ggsave(plot = g_cond_21, filename = paste0(plot_path, "prob_edit_2_given_1.jpg"))

### site 3 given 2
p_3_given_2 = conditional_a_given_b(edit_table = edit_table_by_5, site_number_a = 3, site_number_b = 2)
g_cond_32 = plot_conditional_a_given_b(conditional = p_3_given_2, site_a = "Site3", site_b = "Site2")
g_cond_32
ggsave(plot = g_cond_32, filename = paste0(plot_path, "prob_edit_3_given_2.jpg"))

### site 3 given 1
p_3_given_1 = conditional_a_given_b(edit_table = edit_table_by_5, site_number_a = 3, site_number_b = 1)
g_cond_31 = plot_conditional_a_given_b(conditional = p_3_given_1, site_a = "Site3", site_b = "Site1")
g_cond_31
ggsave(plot = g_cond_31, filename = paste0(plot_path, "prob_edit_3_given_1.jpg"))
