#PURPOSE: Main script for enrichment for candidate genes and pathways

#####Install/load libraries#####
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)
library(ABAEnrichment)
library(ggfortify)
library(gplots)
#####setwd#####
setwd("~/Google\ Drive/BINF\ 4006\ Project\ ")
#####source functions#####
source("~/Google\ Drive/BINF\ 4006\ Project\ /Scripts/Functions.R")
#####pathway dictionary#####
pathlist<-read.delim("~/Google\ Drive/BINF\ 4006\ Project\ /Ensembl2Reactome_All_Levels.txt",sep="\t",stringsAsFactors=F,header=F)
reactome<-pathlist[which(pathlist[,ncol(pathlist)]=="Homo sapiens"),]
reactome_dict<-makeDictionary(reactome,3,1)
path_names_in_dict<-reactome[which(reactome[,3] %in% names(reactome_dict)),4]
#####structure candidates/backgrounds by quantiles#####
struc_cgenes<-lapply(1:414,candidateABAGenes)
struc_bgenes<-lapply(struc_cgenes,FUN=function(x){
  lapply(x,backgroundABAGenes)
  })
#####pathway enrichment#####
candidateGenes<-names(struc_cgenes[[3]][[2]])
backgroundGenes<-struc_bgenes[[3]][[2]]
enrchtest<-NULL
for(i in 1:length(reactome_dict)){
  pathGenes<-reactome_dict[[i]]
  p<-fisher.exact(pathGenes,candidateGenes,backgroundGenes)$p.value
  path_id<-names(reactome_dict)[i]
  path_name<-unique(reactome[which(reactome$V3 %in% path_id),4])
  enrchtest<-rbind(enrchtest,c(path_id,path_name,p))
}
#no significance at all giving random gene sets of equal size. 
#not good
hist(as.numeric(enrchtest[,3]),main="",xlab="") 
#not good-not getting significant pathway enrichments
pi1<-1-pi0est(as.numeric(enrchtest[,3]))$pi0 
#####Permutations#####
i<-1
struc_cgenes<-lapply(1:414,candidateABAGenes)
struc_bgenes<-lapply(struc_cgenes,FUN=function(x){
  lapply(x,backgroundABAGenes)
})
set.seed(1)
#for every pathway, do 1000 permutations with random sets of candidate and
#complementary background genes compute fishers exact test. Check distribution
pathGenes<-reactome_dict[[i]]
candidateGenes<-names(struc_cgenes[[3]][[1]])
i=0
pvals<-c()
while(i<1000){
  candGenes=sample(colnames(mat),length(candidateGenes))
  backgroundGenes=backgroundABAGenes(candidateGenes)
  pvals<-c(pvals,
           fisher.exact(pathGenes,candGenes,backgroundGenes)$p.value)
  i=i+1
}
#####top varying genes#####
vars<-apply(mat,2,var)
sds<-apply(mat,2,sd)
means<-apply(mat,2,mean)
covs<-(sds/means)^2
effectsizes<-apply(mat,2,function(x){(x-means)/sds})
length(sds[which(sds>1.96)])
candidateGenes<-names(sds[which(sds>1.96)])
backgroundGenes<-backgroundABAGenes(candidateGenes)
enrchtest<-NULL
for(i in 1:length(reactome_dict)){
  pathGenes<-reactome_dict[[i]]
  p<-fisher.exact(pathGenes,candidateGenes,backgroundGenes)$p.value
  path_id<-names(reactome_dict)[i]
  path_name<-unique(reactome[which(reactome$V3 %in% path_id),4])
  enrchtest<-rbind(enrchtest,c(path_id,path_name,p))
}
hist(as.numeric(enrchtest[,3]),main="",xlab="") 
#makes total sense
topBrainVaryPaths<-head(enrchtest[order(enrchtest[,3]),2],
     n=length(enrchtest[which(enrchtest[,3]<0.05)]))
#####trying out ABAEnrichment#####
struc_cgenes<-lapply(1:414,candidateABAGenes)
candidateGenes<-rep(1,length(struc_cgenes[[1]][[1]]))
names(candidateGenes)<-names(struc_cgenes[[1]][[1]])
#not good either
res<-aba_enrich(candidateGenes,dataset='adult',cutoff_quantiles =0.01,n_randsets=1000,gene_len = T)

save.image(paste0(getwd(),"/Enrichment.RDa"))
#####pathway x structure matrix#####
#this will yield a spectrum in all brain structures
#for the expression variance of genes in pathways
path_struc_mat<-NULL
path_ids<-c()
for(i in 1:length(reactome_dict)){
  #get genes in pathway
  if(length(reactome_dict[[i]])>5){
    pathGenes<-reactome_dict[[i]]
    path_ids<-c(path_ids,names(reactome_dict)[i])
  }
  else{next}
  #get genes in ABA matrix
  sub_mat<-mat[,which(colnames(mat) %in% pathGenes)]
  #get variance of genes for each structure
  vec<-apply(sub_mat,1,var)
  #in the end should have a pathway by structure matrix
  path_struc_mat<-cbind(path_struc_mat,vec)
}
colnames(path_struc_mat)<-path_ids
#which structures have a high specific path variance?
pathid<-colnames(path_struc_mat)[1]
names<-names(path_struc_mat[order(path_struc_mat[,pathid]),pathid])
inds<-gsub("V","",names)
cat(pathid,"\n")
cat(unique(reactome[which(reactome[,3]==pathid),4]))
head(ontology[inds,6])
#####PCA colored by Pathway expression variation#####
#pick desired pathway(s)
topBrainVaryIds<-unique(reactome[which(reactome[,4] %in% topBrainVaryPaths),3])
for(i in 1:length(topBrainVaryIds)){
  pathwayOfInterest<-topBrainVaryIds[i]
  cat(i," : ",unique(reactome[which(reactome[,3] %in% pathwayOfInterest),4]),"\n")
}
pathwayOfInterestId<-as.character(topEnrichments$aba$path_id[1]) #top enriched pathway in alzheimers
pathwayOfInterestName<-unique(reactome[which(reactome[,3] %in% pathwayOfInterestId),4])
#changing colnames (pathway ids) to something interpretable by ggplot
tmp<-path_struc_mat
strs<-strsplit(colnames(path_struc_mat),"-")
tmp_ids<-unlist(lapply(strs,function(x){paste0(x[1],x[3])}))
rownames(tmp)<-NULL
colnames(tmp)<-tmp_ids
#computing components
pr<-prcomp(mat,
           center = T,
           scale. = T)
#make autoplot params and plot
data<-as.data.frame(tmp)
colorBy<-colnames(data)[which(colnames(path_struc_mat) %in% pathwayOfInterestId)]
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /variability_visualization_pca_example.pdf",width=7,height=4)
autoplot(pr,data=data,
         colour=colorBy,size=1)+
  scale_colour_gradient(name="Pathway variance",
                        low="skyblue",high="red")+
  ggtitle(paste(pathwayOfInterestName))+
  theme(
    title=element_text(face=2,size=10),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=10),
    legend.text=element_text(face=2,size=7)
  )
dev.off()
#####P2 PCA colored by enriched diffeq pathways#####
sgMat<-strucMat(5,"CbCx_Cerebellar Cortex")
spMat<-pathMat(5,"CbCx_Cerebellar Cortex")
pdf(paste0(getwd(),"/Results/Proposal2_PCA_ABAData_CerebellarCortex_ColoredBy1stEnrichedPathway.pdf"),width=6,height=5)
pathwayOfInterestId<-as.character(topEnrichments$aba$path_id[1])
varyingPathwayPCA(pathwayOfInterestId,sgMat,spMat)
dev.off()
pdf(paste0(getwd(),"/Results/Proposal2_PCA_ABAData_CerebellarCortex_ColoredBy2ndEnrichedPathway.pdf"),width=6,height=5)
pathwayOfInterestId<-as.character(topEnrichments$aba$path_id[2])
varyingPathwayPCA(pathwayOfInterestId,sgMat,spMat)
dev.off()
pdf(paste0(getwd(),"/Results/Proposal2_PCA_ABAData_CerebellarCortex_ColoredBy3rdEnrichedPathway.pdf"),width=6,height=5)
pathwayOfInterestId<-as.character(topEnrichments$aba$path_id[3])
varyingPathwayPCA(pathwayOfInterestId,sgMat,spMat)
dev.off()
#####brain structures of interest and their significance#####
load("~/Google Drive/BINF 4006 Project /RDa/DiffExpr.RDa")
#P1
#obs
onts<-NULL
name<-"alz"
data<-path_struc_mat
pheno<-ontology
for(i in 1:nrow(topEnrichments[[name]])){
  pathwayOfInterestId<-as.character(topEnrichments[[name]]$path_id[i]) #top enriched pathway in alzheimers
  pathid<-which(colnames(path_struc_mat) %in% pathwayOfInterestId)
  pathcol<-data[,pathid]
  #ind<-order(pathcol,decreasing=T)
  cut<-quantile(pathcol)["75%"]
  onts[[i]]<-as.data.frame(pheno[pathcol>cut,3:8])
}
strucsWVariablePath<-unlist(lapply(onts,function(x){as.character(x[,2])}))
alz_table<-table(strucsWVariablePath)/length(strucsWVariablePath)
alz_table[order(alz_table,decreasing=F)]*100
sum(alz_table[order(alz_table,decreasing=F)])
#null distributions
alz_names<-names(alz_table)
name<-"alz"
data<-path_struc_mat
pheno<-ontology
alz_null_table<-NULL
for(i in 1:1000){
  nullpaths<-sample(colnames(data),nrow(topEnrichments[[name]]))
  onts<-NULL
  for(j in 1:length(nullpaths)){
    pathwayOfInterestId<-as.character(nullpaths[j]) 
    pathid<-which(colnames(data) %in% pathwayOfInterestId)
    pathcol<-data[,pathid]
    cut<-quantile(pathcol)["75%"]
    onts[[j]]<-as.data.frame(pheno[pathcol>cut,3:8])
  }
  strucsWVariablePath<-unlist(lapply(onts,function(x){as.character(x[,2])}))
  alz_null_table[[i]]<-table(strucsWVariablePath)/length(strucsWVariablePath)
}
#from null table with random structures, get null distribution for null structures
alz_null_struc_freq<-NULL
for(i in 1:length(alz_names)){
  null<-c()
  for(j in 1:length(alz_null_table)){
    null<-c(null,alz_null_table[[j]][names(alz_null_table[[j]]) %in% alz_names[i]])
  }
  alz_null_struc_freq[[i]]<-null
}
names(alz_null_struc_freq)<-alz_names
#compare observed structure proportion with null structure proportion
alz_sig_test_list<-NULL
for(i in 1:length(alz_names)){
  struc<-alz_names[i]
  obs<-alz_table[names(alz_table) %in% struc]
  null_distr<-alz_null_struc_freq[[ i ]]
  alz_sig_test_list[[i]]<-n.test.oneTails(null_distr,obs,0.05)
}
names(alz_sig_test_list)<-alz_names
#plot significance of structures
plot(-log10(unlist(lapply(alz_sig_test_list,function(x){x$p_value}))),ylab="p-values",xlab="structures harboring variability in dz EPs")
abline(h= -log10(0.05),col="red",lwd=2)
#which structure is most significant and what's the frequency of harboring variable pathways?
ind<-which.min(unlist(lapply(alz_sig_test_list,function(x){x$p_value})))
ord_ps<-order(p.adjust(unlist(lapply(alz_sig_test_list,function(x){x$p_value})),method="BH"),decreasing=T)
adj_ps<-p.adjust(unlist(lapply(alz_sig_test_list,function(x){x$p_value})),method="BH")
ord_zs<-order(unlist(lapply(alz_sig_test_list,function(x){x$z_score})),decreasing=T)
alz_sig_test_list[[ind]]
cbind("Z"=lapply(alz_sig_test_list[ord_ps],function(x){x$z_score}),"p-value"=p.adjust(lapply(alz_sig_test_list[ord_ps],function(x){x$p_value}),method="BH"))
cbind("Z"=lapply(alz_sig_test_list[ord_zs],function(x){x$z_score}),"p-value"=adj_ps[ord_zs])
#
#
#P2
onts<-NULL
name<-"pk"
data<-pathMat(5,"CbCx_Cerebellar Cortex")
pheno<-ontology[which(ontology$Structure.Hierarchy.Tier.5=="CbCx_Cerebellar Cortex"),]
for(i in 1:nrow(topEnrichments[[name]])){
  pathwayOfInterestId<-as.character(topEnrichments[[name]]$path_id[i]) #top enriched pathway in alzheimers
  pathid<-which(colnames(data) %in% pathwayOfInterestId)
  pathcol<-data[,pathid]
  #ind<-order(pathcol,decreasing=T)
  cut<-quantile(pathcol)["75%"]
  onts[[i]]<-as.data.frame(pheno[pathcol>cut,3:8])
}
strucsWVariablePath<-unlist(lapply(onts,function(x){as.character(x[,6])}))
pk_table<-table(strucsWVariablePath)/length(strucsWVariablePath)
pk_table[order(pk_table,decreasing=F)]*100
sum(pk_table[order(pk_table,decreasing=F)])

onts<-NULL
name<-"sz"
data<-pathMat(5,"CbCx_Cerebellar Cortex")
pheno<-ontology[which(ontology$Structure.Hierarchy.Tier.5=="CbCx_Cerebellar Cortex"),]
for(i in 1:nrow(topEnrichments[[name]])){
  pathwayOfInterestId<-as.character(topEnrichments[[name]]$path_id[i]) #top enriched pathway in alzheimers
  pathid<-which(colnames(data) %in% pathwayOfInterestId)
  pathcol<-data[,pathid]
  #ind<-order(pathcol,decreasing=T)
  cut<-quantile(pathcol)["75%"]
  onts[[i]]<-as.data.frame(pheno[pathcol>cut,3:8])
}
strucsWVariablePath<-unlist(lapply(onts,function(x){as.character(x[,6])}))
sz_table<-table(strucsWVariablePath)/length(strucsWVariablePath)
sz_table[order(sz_table,decreasing=F)]*100
sum(sz_table[order(sz_table,decreasing=F)])
#
#null distributions
#make null table
pk_names<-names(pk_table)
name<-"pk"
data<-pathMat(5,"CbCx_Cerebellar Cortex")
pheno<-ontology[which(ontology$Structure.Hierarchy.Tier.5=="CbCx_Cerebellar Cortex"),]
pk_null_table<-NULL
for(i in 1:1000){
  nullpaths<-sample(colnames(data),nrow(topEnrichments[[name]]))
  onts<-NULL
  for(j in 1:length(nullpaths)){
    pathwayOfInterestId<-as.character(nullpaths[j]) 
    pathid<-which(colnames(data) %in% pathwayOfInterestId)
    pathcol<-data[,pathid]
    cut<-quantile(pathcol)["75%"]
    onts[[j]]<-as.data.frame(pheno[pathcol>cut,3:8])
  }
  strucsWVariablePath<-unlist(lapply(onts,function(x){as.character(x[,6])}))
  pk_null_table[[i]]<-table(strucsWVariablePath)/length(strucsWVariablePath)
}
#from null table with random structures, get null distribution for null structures
pk_null_struc_freq<-NULL
for(i in 1:length(pk_names)){
  null<-c()
  for(j in 1:length(pk_null_table)){
    null<-c(null,pk_null_table[[j]][names(pk_null_table[[j]]) %in% pk_names[i]])
  }
  pk_null_struc_freq[[i]]<-null
}
names(pk_null_struc_freq)<-pk_names
#compare observed structure proportion with null structure proportion
pk_sig_test_list<-NULL
for(i in 1:length(pk_names)){
  struc<-pk_names[i]
  obs<-pk_table[names(pk_table) %in% struc]
  null_distr<-pk_null_struc_freq[[ which(names(pk_null_struc_freq) %in% struc) ]]
  pk_sig_test_list[[i]]<-n.test.oneTails(null_distr,obs,0.05)
}
names(pk_sig_test_list)<-pk_names
#plot significance of structures
plot(-log10(unlist(lapply(pk_sig_test_list,function(x){x$p_value}))),ylab="p-values",xlab="structures harboring variability in dz EPs")
abline(h= -log10(0.05),col="red",lwd=2)
#which structure is most significant and what's the frequency of harboring variable pathways?
ind<-which.min(unlist(lapply(pk_sig_test_list,function(x){x$p_value})))
ord_ps<-order(p.adjust(unlist(lapply(pk_sig_test_list,function(x){x$p_value})),method="BH"),decreasing=T)
adj_ps<-p.adjust(unlist(lapply(pk_sig_test_list,function(x){x$p_value})),method="BH")
ord_zs<-order(unlist(lapply(pk_sig_test_list,function(x){x$z_score})),decreasing=T)
pk_sig_test_list[[ind]]
cbind("Z"=lapply(pk_sig_test_list[ord_ps],function(x){x$z_score}),"p-value"=p.adjust(lapply(pk_sig_test_list[ord_ps],function(x){x$p_value}),method="BH"))
cbind("Z"=lapply(pk_sig_test_list[ord_zs],function(x){x$z_score}),"p-value"=adj_ps[ord_zs])
#
#
sz_names<-names(sz_table)
name<-"sz"
data<-pathMat(5,"CbCx_Cerebellar Cortex")
pheno<-ontology[which(ontology$Structure.Hierarchy.Tier.5=="CbCx_Cerebellar Cortex"),]
sz_null_table<-NULL
for(i in 1:1000){
  nullpaths<-sample(colnames(data),nrow(topEnrichments[[name]]))
  onts<-NULL
  for(j in 1:length(nullpaths)){
    pathwayOfInterestId<-as.character(nullpaths[j]) 
    pathid<-which(colnames(data) %in% pathwayOfInterestId)
    pathcol<-data[,pathid]
    cut<-quantile(pathcol)["75%"]
    onts[[j]]<-as.data.frame(pheno[pathcol>cut,3:8])
  }
  strucsWVariablePath<-unlist(lapply(onts,function(x){as.character(x[,6])}))
  sz_null_table[[i]]<-table(strucsWVariablePath)/length(strucsWVariablePath)
}
#from null table with random structures, get null distribution for null structures
sz_null_struc_freq<-NULL
for(i in 1:length(sz_names)){
  null<-c()
  for(j in 1:length(sz_null_table)){
    null<-c(null,sz_null_table[[j]][names(sz_null_table[[j]]) %in% sz_names[i]])
  }
  sz_null_struc_freq[[i]]<-null
}
names(sz_null_struc_freq)<-sz_names
#compare observed structure proportion with null structure proportion
sz_sig_test_list<-NULL
for(i in 1:length(sz_names)){
  struc<-sz_names[i]
  obs<-sz_table[names(sz_table) %in% struc]
  null_distr<-sz_null_struc_freq[[ which(names(sz_null_struc_freq) %in% struc) ]]
  sz_sig_test_list[[i]]<-n.test.oneTails(null_distr,obs,0.05)
}
names(sz_sig_test_list)<-sz_names
#plot significance of structures
plot(-log10(unlist(lapply(sz_sig_test_list,function(x){x$p_value}))),ylab="p-values",xlab="structures harboring variability in dz EPs")
abline(h= -log10(0.05),col="red",lwd=2)
#which structure is most significant and what's the frequency of harboring variable pathways?
ind<-which.min(unlist(lapply(sz_sig_test_list,function(x){x$p_value})))
ord_ps<-order(p.adjust(unlist(lapply(sz_sig_test_list,function(x){x$p_value})),method="BH"),decreasing=T)
adj_ps<-p.adjust(unlist(lapply(sz_sig_test_list,function(x){x$p_value})),method="BH")
ord_zs<-order(unlist(lapply(sz_sig_test_list,function(x){x$z_score})),decreasing=T)
sz_sig_test_list[[ind]]
cbind("Z"=lapply(sz_sig_test_list[ord_ps],function(x){x$z_score}),"p-value"=p.adjust(lapply(sz_sig_test_list[ord_ps],function(x){x$p_value}),method="BH"))
cbind("Z"=lapply(sz_sig_test_list[ord_zs],function(x){x$z_score}),"p-value"=adj_ps[ord_zs])
#####trying out GSEABase#####
source("https://bioconductor.org/biocLite.R")
biocLite("GSEABase")
library(GSEABase)
#####save image#####
save.image("~/Documents/Columbia/Courses/TRANSLATIONAL_BIOINFORMATICS/project/RDa/Enrichment.RDa")
