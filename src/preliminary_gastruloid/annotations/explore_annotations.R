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

setwd("~/typewriter_analysis")


annotations <- read.csv("data/mGASv2_Lane2_Group1_cell_annotation.csv")
annotations <- annotations[,2:(ncol(annotations))]
annotations <- data.frame(annotations)

plot1 <- ggplot(data=annotations, aes(x=jax_major_trajectory)) + geom_bar() + xlab("Jax Annotation") + ylab("Total counts") +  theme(legend.title = element_blank(),axis.text.x=element_text(angle = 90))
ggsave("jax_annotations.pdf",path="results/preliminary_gastruloid/annotations/", width=25,height= 18, units = "cm")


plot2 <- ggplot(data=annotations, aes(x=pijuan_celltype)) + geom_bar() + xlab("Pijuan Annotation") + ylab("Total counts") +  theme(legend.title = element_blank(),axis.text.x=element_text(angle = 90))
ggsave("pijuan_annotations.pdf",path="results/preliminary_gastruloid/annotations/", width=25,height= 18, units = "cm")

#create a manual annotation to match the hierarchy in this preprint:
#https://www.biorxiv.org/content/10.1101/2022.11.01.514697v1.full.pdf 
# based on pijaun annotations: we merge pijuan Brain annotations with SC
annotations_xml <- c()
for(i in 1:length(annotations$pijuan_celltype)) {
  if(annotations$pijuan_celltype[i] == "Paraxial mesoderm") {
    annotations_xml <- c(annotations_xml,"PxMD")
  }
  else if(annotations$pijuan_celltype[i] == "Pharyngeal mesoderm") {
    annotations_xml <- c(annotations_xml,"PhMD")
  }
  else if(annotations$pijuan_celltype[i] == "Somitic mesoderm") {
    annotations_xml <- c(annotations_xml,"SMD")
  }
  else if(annotations$pijuan_celltype[i] == "Spinal chord") {
    annotations_xml <- c(annotations_xml,"SC")
  }
  else if(annotations$pijuan_celltype[i] == "Forebrain/Midbrain/Hindbrain") {
    annotations_xml <- c(annotations_xml,"SC")
  }
  else if(annotations$pijuan_celltype[i] == "NMP") {
    annotations_xml <- c(annotations_xml,"NMP")
  }
  else {annotations_xml <- c(annotations_xml,NA)}
}

length(annotations_xml)
sum(is.na(annotations_xml))
#iterate through the pijuan annotations


taxa <- paste0(annotations[which(! is.na(annotations_xml)),]$Cell,"-1")
data <- annotations_xml[which(! is.na(annotations_xml))]
data_taxa <- data.frame(taxa,data)

write.table(data_taxa,file="results/preliminary_gastruloid/annotations/annotations.csv",sep=", ",dec = " ",quote=F,row.names=F,col.names=F)
