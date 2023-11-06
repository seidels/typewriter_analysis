## ---------------------------
##
## Script name: explore_annotations
##
## Purpose of script: explore cell-type annotations for the gastruloid analysis and create annotations file for the BDMM
##
## Author: Antoine Zwaans
##
## Date Created: 2023-09-01
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##

require(tidyverse)
require(data.table)

setwd("~/Projects/typewriter_analysis")

#input and filter datasets by set of 8 targetBCs
edit_table_by_6 <- read.csv("data/preliminary_gastruloid/mGASv2_Lane2_CellByTape_10X_bamExtractV2_t3_collapse.csv", stringsAsFactors = F, header = T,na.strings=c("","NA"))
edit_table_by_5 <- edit_table_by_6[,1:7]

filters <- read.csv("data/preliminary_gastruloid/mGASv2_TargetBC_selected8.csv", header = T,na.strings=c("","NA"))

filtered_dataset <- c()
for(i in filters$TargetBC) {
  filtered_dataset <- rbind(filtered_dataset,edit_table_by_5[which(edit_table_by_5$TargetBC == i),])
}

filtered_by_targetBC <- filtered_dataset

frequent_target_bcs <- filters$TargetBC
##further filter by cell that have ALL 8. 
full_filter<- filtered_dataset
for (cell in unique(full_filter$Cell)){
  
  cell_target_bcs = full_filter[which(full_filter$Cell == cell), "TargetBC"]
  keep_cell = all(is.element(el = frequent_target_bcs, set = cell_target_bcs))
  
  if (! keep_cell){
    full_filter = full_filter[which(!(full_filter$Cell == cell)), ]
  }
}



#input annotations
annotations <- read.csv("data/mGASv2_Lane2_Group1_cell_annotation.csv")
annotations <- annotations[,2:(ncol(annotations))]

annotations <- data.frame(annotations)

#checking that the cells in this annotation table correspond to the 780 subset that 
#choi and sam sent us. filter out other cells

annotations$Cell <- paste0(annotations$Cell, "-1") 

cell_names_filtered780 <- full_filter$Cell

cells_to_keep <- Reduce(intersect,list(annotations$Cell,cell_names_filtered780))

annotations_filtered <- annotations[annotations$Cell %in% cells_to_keep,]

nrow(annotations_filtered)
#there are only 296 cells in the annotatons file that match the 780 

plot1 <- ggplot(data=annotations_filtered, aes(x=jax_major_trajectory)) + geom_bar() + xlab("Jax Annotation") + ylab("Total counts") +  theme(legend.title = element_blank(),axis.text.x=element_text(angle = 90))
ggsave("jax_annotations.pdf",path="results/preliminary_gastruloid/annotations/", width=25,height= 18, units = "cm")


plot2 <- ggplot(data=annotations_filtered, aes(x=pijuan_celltype)) + geom_bar() + xlab("Pijuan Annotation") + ylab("Total counts") +  theme(legend.title = element_blank(),axis.text.x=element_text(angle = 90))
ggsave("pijuan_annotations.pdf",path="results/preliminary_gastruloid/annotations/", width=25,height= 18, units = "cm")

#create a manual annotation to match the hierarchy in this preprint:
#https://www.biorxiv.org/content/10.1101/2022.11.01.514697v1.full.pdf
# based on pijaun annotations: we merge pijuan Brain annotations with SC

length(unique(annotations_filtered$pijuan_celltype))
annotations_xml <- c()
for(i in 1:length(annotations_filtered$pijuan_celltype)) {
  if(annotations_filtered$pijuan_celltype[i] == "Paraxial mesoderm") {
    annotations_xml <- c(annotations_xml,"PxMD")
  }
  else if(annotations_filtered$pijuan_celltype[i] == "Pharyngeal mesoderm") {
    annotations_xml <- c(annotations_xml,"PhMD")
  }
  else if(annotations_filtered$pijuan_celltype[i] == "Somitic mesoderm") {
    annotations_xml <- c(annotations_xml,"SMD")
  }
  else if(annotations_filtered$pijuan_celltype[i] == "Spinal chord") {
    annotations_xml <- c(annotations_xml,"SC")
  }
  else if(annotations_filtered$pijuan_celltype[i] == "Forebrain/Midbrain/Hindbrain") {
    annotations_xml <- c(annotations_xml,"SC")
  }
  else if(annotations_filtered$pijuan_celltype[i] == "NMP") {
    annotations_xml <- c(annotations_xml,"NMP")
  }
  else {annotations_xml <- c(annotations_xml,NA)}
}

length(annotations_xml)
sum(is.na(annotations_xml))

unique(annotations_xml)
#iterate through the pijuan annotations


taxa <- annotations_filtered$Cell[which(! is.na(annotations_xml))]
data <- annotations_xml[which(! is.na(annotations_xml))]
data_taxa <- data.frame(taxa,data)

write.table(data_taxa,file="results/preliminary_gastruloid/annotations/annotations.csv",sep=", ",dec = " ",quote=F,row.names=F,col.names=F)
