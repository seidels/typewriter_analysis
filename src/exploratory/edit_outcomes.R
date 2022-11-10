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

setwd("~/Projects/typewriter_analysis/")      # Sophie's working directory (mac)
library(dplyr)
library(ggplot2)
library(tidyr)
library(boot)
library(cowplot)


## ---------------------------

# load the preprocessed data ---------

edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))

# basic properties of the dataset ----


# total number of cells: 16863 doesn't match the 12000 mentioned in paper?

length(unique(edit_table_by_5$Cell))

#total number of unique TargetBC: 17 + NA = 18

length(unique(edit_table_by_5$TargetBC))

#total unique edit sequences at all positions included: 19 + NA:

length(unique(unlist(edit_table_by_5[,3:7])))

#total cells:



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

#testing whether insert distributions are independent of sites?

chisq.test(t(table(insert_count_site)))

ggplot(insert_count_site, aes(x=site,fill=sequence)) + geom_bar() + xlab("Site") + ylab("Total counts")+ theme(legend.title = element_blank())
ggsave("bulk_insert_count_site_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

#same figure without coloured insert ditributions
ggplot(insert_count_site, aes(x=site)) + geom_bar() + xlab("Site") + ylab("Total counts")+ theme(legend.title = element_blank())
ggsave("bulk_count_site_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

counts_per_site <- count(insert_count_site,site)

diffs <- c()
for(i in 2:length(counts_per_site$n)-1) {
  diffs <- c(diffs,(counts_per_site$n[i+1]-counts_per_site$n[i])/counts_per_site$n[i])

}
diffs <- (diffs* (-1)) *100

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
insert_count_targBC <- data.frame(edit_table_by_5[,c(1,2,3,4,5,6,7)]) %>% gather(targetBC, sequence,3:7)
insert_count_targBC <- insert_count_targBC[- which(is.na(insert_count_targBC$sequence)),]


#testing whether insert distributions are independent of sites?
chisq.test(t(table(insert_count_targBC)))


insert_count_targBC$sequence <- substring(insert_count_targBC$sequence, 1,3)
ggplot(insert_count_targBC, aes(x=TargetBC,fill=sequence)) + geom_bar() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("TargetBC") + ylab("Count") + theme(legend.title = element_blank())
ggsave("insert_count_targetBC_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


cell_count_targBC <- count(insert_count_targBC, insert_count_targBC$TargetBC,insert_count_targBC$Cell)
unique_cells <- c()
for( i in unique(cell_count_targBC$`insert_count_targBC$TargetBC`)) {
  unique_cells <- c(unique_cells,length(unique(cell_count_targBC[cell_count_targBC$`insert_count_targBC$TargetBC` == i,2])))

}

cells_per_targetBC <- data.frame(targetBC= unique(cell_count_targBC$`insert_count_targBC$TargetBC`),unique_cells = unique_cells )
ggplot(cells_per_targetBC,aes(x=targetBC,y=unique_cells)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) +
  xlab("TargetBC") + ylab("Number of cells")
ggsave("cells_per_targetBC.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")
## divide by total number of integrations (not cells)


# insert frequencies - TargetBC

insert_freq_targBC <- data.frame(edit_table_by_5[,c(1,2,3,4,5,6,7)]) %>% gather(site, sequence,3:7)
insert_freq_targBC$sequence <- substring(insert_freq_targBC$sequence, 1,3)
insert_freq_targBC <- count(insert_freq_targBC, insert_freq_targBC$TargetBC, insert_freq_targBC$sequence)
colnames(insert_freq_targBC) <- c("target_BC","sequence","count")
ggplot(insert_freq_targBC, aes(fill=sequence, y=count, x=target_BC)) + theme(axis.text.x = element_text(angle = 90)) +
  geom_bar(position="fill", stat="identity") +  theme(legend.title = element_blank()) +
  xlab("TargetBC") + ylab("Frequency")  + theme(legend.title = element_blank())
ggsave("insert_freq_targetBC.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

# insert frequencies - TargetBC - no NA

insert_freq_targBC <- data.frame(edit_table_by_5[,c(2,3,4,5,6,7)]) %>% gather(targetBC, sequence,2:6)
insert_freq_targBC <- insert_freq_targBC[ which(! is.na(insert_freq_targBC$sequence)),]
insert_freq_targBC$sequence <- substring(insert_freq_targBC$sequence, 1,3)
insert_freq_targBC <- count(insert_freq_targBC, insert_freq_targBC$TargetBC, insert_freq_targBC$sequence)
colnames(insert_freq_targBC) <- c("target_BC","sequence","count")
ggplot(insert_freq_targBC, aes(fill=sequence, y=count, x=target_BC)) +
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90)) +
  xlab("TargetBC") + ylab("Frequency") + theme(legend.title = element_blank())
ggsave("insert_freq_targetBC_no_na.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")


insert_freq_targBC_cell <- data.frame(edit_table_by_5[,c(1,2,3,4,5,6,7)]) %>% gather(TargetBC, site,3:7)


## ---------------------------

## Check edit outcome probability conditional on previous site

edit_table_by_5 = read.csv("data/Supplementary_File_2_DataTableMOI19.csv", stringsAsFactors = F, header = T, na.strings=c("","NA"))
edit_table_by_5[is.na(edit_table_by_5)] <- "NA"
all_inserts <- unique(unlist(edit_table_by_5[,3:7]))
contiguous_inserts <- matrix(0,nrow = length(all_inserts),ncol = length(all_inserts))
rownames(contiguous_inserts) = all_inserts
colnames(contiguous_inserts) = all_inserts

#counting manually the joint probs.
for(i in 1:(length(edit_table_by_5$Cell))) {
  for (j in 3:6) {
    for (a in 1:length(all_inserts)) {
      for (b in 1:length(all_inserts)) {
        if((edit_table_by_5[i,j] == all_inserts[a]) && (edit_table_by_5[i,j+1] == all_inserts[b])) {
          contiguous_inserts[a,b] = contiguous_inserts[a,b] + 1
        }
      }
    }
  }



}

#converting to conditional probs
bigram_matrix <- data.frame(contiguous_inserts)
bigram_matrix <- bigram_matrix[-which(rownames(bigram_matrix) == "NA"), -which(colnames(bigram_matrix) == "NA.")]
#divide rows by their rowsums:
bigram_matrix <- apply(bigram_matrix, 1, function(x) {x=x/sum(x)})

#very dirty way to melt the matrix to a 3 column thing
count <- c()
position1 <- c()
for( i in rownames(bigram_matrix)) {

  position1 <- c(position1,rep(i,19))

}
position1 <- unlist(position1)
position2 <- rep(colnames(bigram_matrix),19)
counts <- as.vector(bigram_matrix)

#clean up insert symbols
position1 <- substring(position1, 1,3)
position2 <- substring(position2, 1,3)

heatmap_data <- data.frame(pos1 = position1, pos2=position2, count = counts)
View(heatmap_data)
ggplot(heatmap_data, aes(x=pos1, y=pos2, fill=count)) +
  geom_tile() + xlab("Position 1") + ylab("Position 2") +
  theme(legend.title = element_text("Normalized Count"))
+ geom_tile(color = "black") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") + labs(fill="P(Position2|Position1)") + coord_fixed()

ggsave("conditional_probability_matrix.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

## bootstrap

# Bootstrap 95% CI for the difference in edit counts per site

differences <- function(data, indices) {
  d <- data[indices,]
  insert_count_site <- data.frame(d[,3:7]) %>% gather(site, sequence,1:5)
  insert_count_site$sequence <- substring(insert_count_site$sequence, 1,3)
  insert_count_site <- insert_count_site[- which(is.na(insert_count_site$sequence)),]
  counts <- count(insert_count_site,site)
  count_diff <- c()
  for (i in 1:length(counts$n)-1) {
    count_diff <- c(count_diff,(counts$n[i]-counts$n[i+1])/counts$n[i])
  }
  return(count_diff)
}

results <- boot(data=edit_table_by_5, statistic=differences,R=1000)
steps_bootstrap <- data.frame(results$t)
colnames(steps_bootstrap) <- c("step_1_2","step_2_3","step_3_4","step_4_5")
steps_bootstrap <- gather(steps_bootstrap,Step)
p0 <- ggplot(data=steps_bootstrap,aes(x=value*100,fill=Step)) + geom_histogram(aes(y=..count../sum(..count..)),bins=250) + xlab("Percent change edited sites") + ylab("Density")+ theme(legend.position = "none")
ggsave("bootstrap_step_norm.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

data_1_2 <- steps_bootstrap[steps_bootstrap$Step=="step_1_2",]
p1 <- ggplot(data=data_1_2,aes(x=value*100)) + geom_histogram(aes(y=..count../sum(..count..)),bins=100,color="coral") + xlab("Percent change sites 1 to 2") + ylab("Density")+
theme(legend.position = "none") + geom_vline(xintercept = mean(data_1_2$value)*100) +
geom_vline(xintercept = sort(data_1_2$value)[975]*100,linetype = "dashed") +
geom_vline(xintercept = sort(data_1_2$value)[25]*100,linetype = "dashed")

data_2_3 <- steps_bootstrap[steps_bootstrap$Step=="step_2_3",]
p2 <- ggplot(data=data_2_3,aes(x=value*100)) + geom_histogram(aes(y=..count../sum(..count..)),bins=100,color="limegreen") + xlab("Percent change sites 2 to 3") + ylab("Density")+
  theme(legend.position = "none") + geom_vline(xintercept = mean(data_2_3$value)*100) +
  geom_vline(xintercept = sort(data_2_3$value)[975]*100,linetype = "dashed") +
  geom_vline(xintercept = sort(data_2_3$value)[25]*100,linetype = "dashed")

data_3_4 <- steps_bootstrap[steps_bootstrap$Step=="step_3_4",]
p3 <- ggplot(data=data_3_4,aes(x=value*100)) + geom_histogram(aes(y=..count../sum(..count..)),bins=100,color="turquoise") + xlab("Percent change sites 3 to 4") + ylab("Density")+
  theme(legend.position = "none") + geom_vline(xintercept = mean(data_3_4$value)*100) +
  geom_vline(xintercept = sort(data_3_4$value)[975]*100,linetype = "dashed") +
  geom_vline(xintercept = sort(data_3_4$value)[25]*100,linetype = "dashed")

data_4_5 <- steps_bootstrap[steps_bootstrap$Step=="step_4_5",]
p4 <- ggplot(data=data_4_5,aes(x=value*100)) + geom_histogram(aes(y=..count../sum(..count..)),bins=100,color="magenta") + xlab("Percent change sites 4 to 5") + ylab("Density")+
  theme(legend.position = "none") + geom_vline(xintercept = mean(data_4_5$value)*100) +
  geom_vline(xintercept = sort(data_4_5$value)[975]*100,linetype = "dashed") +
  geom_vline(xintercept = sort(data_4_5$value)[25]*100,linetype = "dashed")


bottom_row <- plot_grid(p1, p2, p3,p4, nrow=1)
plot_grid(p0, bottom_row, nrow=2)
ggsave("bootstrap_detailed.pdf",path="/Users/azwaans/typewriter_analysis/results/exploratory", width=25,height= 18, units = "cm")

