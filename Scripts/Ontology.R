#PURPOSE: T.B. project-compute PCA for data

#####Install/load libraries#####
source("https://bioconductor.org/biocLite.R")
biocLite("ABAEnrichment")
library(ABAEnrichment)
#####load Data#####
load("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/ABAprocessing.RDa")
#####Ontology#####
struc_ids = paste("Allen:", unique(dataset_adult$structure), sep="")
names<-get_name(struc_ids)
superstrucs<-lapply(names(names),get_superstructures)
superstruc_names<-lapply(superstrucs,get_name)
Superstructures1<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][1]}))
Superstructures2<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][2]}))
Superstructures3<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][3]}))
Superstructures4<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][4]}))
Superstructures5<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][5]}))
Superstructures6<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][6]}))
Superstructures7<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][7]}))
Superstructures8<-unlist(lapply(1:length(superstruc_names),function(x){superstruc_names[[x]][8]}))
ontology<-data.frame("Structure.Hierarchy.Tier.1"=Superstructures1,"Structure.Hierarchy.Tier.2"=Superstructures2,"Structure.Hierarchy.Tier.3"=Superstructures3,"Structure.Hierarchy.Tier.4"=Superstructures4,"Structure.Hierarchy.Tier.5"=Superstructures5,"Structure.Hierarchy.Tier.6"=Superstructures6,"Structure.Hierarchy.Tier.7"=Superstructures7,"Structure.Hierarchy.Tier.8"=Superstructures8)
#####save image#####
save.image("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/Ontology.RDa")
