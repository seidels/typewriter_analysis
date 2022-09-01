## ---------------------------
##
## Script name: edit_outcomes
##
## Purpose of script: Analyze the editing outcomes of the lineage tracing experiment in https://doi.org/10.1038/s41586-022-04922-8
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

setwd("/Users/azwaans/typewriter_analysis/")      # Sophie's working directory (mac)
setwd("~/Projects/typewriter_analysis/")

## ---------------------------

## libraries
library(ggplot2)
library(reshape2)

## ---------------------------

## figure info
plot_path = "results/exploratory/"

# load the preprocessed data ---------

edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))

# basic properties of the dataset ----
# total number of cells: 16863

length(unique(edit_table_by_5$Cell))

#total number of unique TargetBC: 17 + NA = 18

length(unique(edit_table_by_5$TargetBC))

#total unique edit sequences at all positions included: 19 + NA:

length(unique(unlist(edit_table_by_5[,3:7])))


## Check edit outcome probability site-independent


## histograms of insert frequencies/counts -------------------------

#insert counts - bulk
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]),useNA = "always")) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$Var1[1] <- "NA"
bulk_insert_count$Var1 <- factor(bulk_insert_count$Var1,levels=bulk_insert_count$Var1)
ggplot(data=bulk_insert_count, aes(x=Var1,y=Freq,fill=Var1)) + geom_col() + xlab("Insert sequence") + ylab("Total counts") +  theme(legend.title = element_blank())
ggsave("bulk_insert_count.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

#insert counts - bulk - no NA
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]))) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$Var1 <- factor(bulk_insert_count$Var1,levels=bulk_insert_count$Var1)
ggplot(data=bulk_insert_count, aes(x=Var1,y=Freq,fill=Var1)) + geom_col() + xlab("Insert sequence") + ylab("Total counts") + theme(legend.title = element_blank())
ggsave("bulk_insert_count_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

#insert frequencies - bulk
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]),useNA = "always")) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$Var1[1] <- "NA"
bulk_insert_count$Var1 <- factor(bulk_insert_count$Var1,levels=bulk_insert_count$Var1)
ggplot(data=bulk_insert_count, aes(x=Var1,y=Freq/sum(Freq),fill=Var1)) + geom_col() + xlab("Insert sequence") + ylab("Frequency")  +  theme(legend.title = element_blank())
ggsave("bulk_insert_freq.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

#insert frequencies - bulk  - no NA
bulk_insert_count <- data.frame(table(unlist(edit_table_by_5[,3:7]))) %>% arrange(desc(Freq))
bulk_insert_count$Var1 <- substring(bulk_insert_count$Var1, 1,3)
bulk_insert_count$Var1 <- factor(bulk_insert_count$Var1,levels=bulk_insert_count$Var1)
ggplot(data=bulk_insert_count, aes(x=Var1,y=Freq/(sum(Freq)),fill=Var1)) + geom_col() + xlab("Insert sequence") + ylab("Total counts") +  theme(legend.title = element_blank())
ggsave("bulk_insert_freq_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


#insert counts - site
insert_count_site <- data.frame(edit_table_by_5[,3:7]) %>% gather(site, sequence,1:5)
insert_count_site$sequence <- substring(insert_count_site$sequence, 1,3)
ggplot(insert_count_site, aes(x=site,fill=sequence)) + geom_histogram(stat="count") + xlab("Site") + ylab("Total counts") + theme(legend.title = element_blank())
ggsave("bulk_insert_count_site.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


#insert counts - site - no NA
insert_count_site <- data.frame(edit_table_by_5[,3:7]) %>% gather(site, sequence,1:5)
insert_count_site$sequence <- substring(insert_count_site$sequence, 1,3)
insert_count_site <- insert_count_site[- which(is.na(insert_count_site$sequence)),]
ggplot(insert_count_site, aes(x=site,fill=sequence)) + geom_bar() + xlab("Site") + ylab("Total counts")+ theme(legend.title = element_blank())
ggsave("bulk_insert_count_site_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


# insert frequency - site
insert_freq_site  <- data.frame(edit_table_by_5[,3:7]) %>% gather(site, sequence,1:5)
insert_freq_site$sequence <- substring(insert_freq_site$sequence, 1,3)
insert_freq_site <- count(insert_freq_site, insert_freq_site$site, insert_freq_site$sequence)
colnames(insert_freq_site) <- c("site","sequence","count")
ggplot(insert_freq_site, aes(fill=sequence, y=count, x=site)) +
  geom_bar(position="fill", stat="identity") + xlab("Site") + ylab("Frequency") + theme(legend.title = element_blank())
ggsave("bulk_insert_freq_site.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")



# insert frequency - site - no na
insert_freq_site  <- data.frame(edit_table_by_5[,3:7]) %>% gather(site, sequence,1:5)
insert_freq_site$sequence <- substring(insert_freq_site$sequence, 1,3)
insert_freq_site <- insert_freq_site[- which(is.na(insert_freq_site$sequence)),]
insert_freq_site <- count(insert_freq_site, insert_freq_site$site, insert_freq_site$sequence)
colnames(insert_freq_site) <- c("site","sequence","count")
ggplot(insert_freq_site, aes(fill=sequence, y=count, x=site)) +
  geom_bar(position="fill", stat="identity") + xlab("Site") + ylab("Frequency") + theme(legend.title = element_blank())
ggsave("bulk_insert_freq_site_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

# insert counts - TargetBC
insert_count_targBC <- data.frame(edit_table_by_5[,c(2,3,4,5,6,7)]) %>% gather(targetBC, sequence,2:6)
insert_count_targBC <- insert_count_targBC[, c(1,3)]
insert_count_targBC$sequence <- substring(insert_count_targBC$sequence, 1,3)
ggplot(insert_count_targBC, aes(x=TargetBC,fill=sequence)) + geom_bar() + theme(axis.text.x = element_text(angle = 90)) +
xlab("TargetBC") + ylab("Count") + theme(legend.title = element_blank())
ggsave("insert_count_targetBC.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

# insert counts - TargetBC - no na
insert_count_targBC <- data.frame(edit_table_by_5[,c(2,3,4,5,6,7)]) %>% gather(targetBC, sequence,2:6)
insert_count_targBC <- insert_count_targBC[, c(1,3)]
insert_count_targBC <- insert_count_targBC[- which(is.na(insert_count_targBC$sequence)),]
insert_count_targBC$sequence <- substring(insert_count_targBC$sequence, 1,3)
ggplot(insert_count_targBC, aes(x=TargetBC,fill=sequence)) + geom_bar() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("TargetBC") + ylab("Count") + theme(legend.title = element_blank())
ggsave("insert_count_targetBC_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


# insert frequencies - TargetBC
insert_freq_targBC <- data.frame(edit_table_by_5[,c(2,3,4,5,6,7)]) %>% gather(targetBC, sequence,2:6)
insert_freq_targBC <- insert_freq_targBC[, c(1,3)]
insert_freq_targBC$sequence <- substring(insert_freq_targBC$sequence, 1,3)
insert_freq_targBC <- count(insert_freq_targBC, insert_freq_targBC$TargetBC, insert_freq_targBC$sequence)
colnames(insert_freq_targBC) <- c("target_BC","sequence","count")
ggplot(insert_freq_targBC, aes(fill=sequence, y=count, x=target_BC)) + theme(axis.text.x = element_text(angle = 90)) +
geom_bar(position="fill", stat="identity") +  theme(legend.title = element_blank()) +
xlab("TargetBC") + ylab("Frequency")  + theme(legend.title = element_blank())
ggsave("insert_freq_targetBC.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


# insert frequencies - TargetBC - no NA
insert_freq_targBC <- data.frame(edit_table_by_5[,c(2,3,4,5,6,7)]) %>% gather(targetBC, sequence,2:6)
insert_freq_targBC <- insert_freq_targBC[, c(1,3)]
insert_freq_targBC <- insert_freq_targBC[ which(! is.na(insert_freq_targBC$sequence)),]
insert_freq_targBC$sequence <- substring(insert_freq_targBC$sequence, 1,3)
insert_freq_targBC <- count(insert_freq_targBC, insert_freq_targBC$TargetBC, insert_freq_targBC$sequence)
colnames(insert_freq_targBC) <- c("target_BC","sequence","count")
ggplot(insert_freq_targBC, aes(fill=sequence, y=count, x=target_BC)) +
geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90)) +
xlab("TargetBC") + ylab("Frequency") + theme(legend.title = element_blank())
ggsave("insert_freq_targetBC_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


## ---------------------------
## explore edit outcomes

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
