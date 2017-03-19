#PURPOSE: T.B. project-Analyze, curate and anayze alzheimers geo data

#####Download and load packages#####
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("ABAData")
biocLite("ABAEnrichment")
biocLite("GEOquery")
install.packages("data.table")
biocLite("limma")
install.packages("ggfortify")
install.packages("rgl")
biocLite("sva")
#ABAData has been out for like 13 months but ABAEnrichment has only been out for <1 year.
require(ABAData)
require(ABAEnrichment)
require(GEOquery)
require(data.table)
require(ggfortify)
require(limma)
require(sva)
require(biomaRt)
#####Alzheimers data with sva and pca with ABAData#####
#get datasets and pheno data
gsem=getGEO("GSE48350",GSEMatrix = T)
pheno<-pData(phenoData(gsem[[1]]))
gsenom= getGEO("GSE48350",GSEMatrix = F);
names(GPLList(gsenom))
#potential batch effect contributing variables:
names(which(apply(pheno,2,function(x){length(table(x))>1})))
#making dataset
probesets <- Table(GPLList(gsenom)[[1]])$ID
data.matrix <- do.call('cbind',lapply(GSMList(gsenom),function(x){
  tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
#get sample batch variables
agetmp<-unlist(lapply(strsplit(as.character(pheno$title),"_"),function(x){x[3]}))
age<-gsub("yrs","",agetmp)
gender<-ifelse(grepl("_female",pheno$title),"female","male")
regions<-unlist(lapply(strsplit(as.character(pheno$title),"_"),function(x){x[1]}))
proc_regions<-gsub("-","",gsub(" ","",tolower(regions)))
dz<-ifelse(grepl("_AD_",pheno$title),"AD","ctrl")
stage<-unlist(lapply(strsplit(as.character(pheno$characteristics_ch1.3)," "),function(x){x[3]}))
apoe<-unlist(lapply(strsplit(as.character(pheno$characteristics_ch1.4)," "),function(x){x[3]}))
alzpheno<-data.frame(
  gender,
  proc_regions,
  dz,
  age,
  stage,
  apoe
)
#try to get this to work to correct for "unknown and unwanted" batches
ind<-round(runif(round(ncol(data.matrix)*0.1,0),1,ncol(data.matrix)),0)
train<-data.matrix[,ind]
np<-alzpheno[ind,]
test<-data.matrix[,-c(ind)]
#making model matrices for batch correction with sva
mod<-model.matrix(~dz+proc_regions, data=np)
mod0 = model.matrix(~gender,data=np) 
#may need to add adjustment variables
#sva removes batch effects and unewanted variation-but is primarily
#used for D.E.
#S.V.s are covariates from HD data
#used if you don't know the batches
n.sv = num.sv(train,mod,method="leek") #why 45?
svobj = sva(train,mod,mod0,n.sv=n.sv)
#computing surrogates based on training data then predicting
#G.E. levels for testing data
fsv1<-fsva(train,mod,svobj,test,method="exact")
fsv2<-fsva(train,mod,svobj,data.matrix,method="exact")
#if you know the batches
mod<-model.matrix(~dz+proc_regions, data=alzpheno)
combat_data = ComBat(dat=data.matrix, 
                     batch=as.integer(as.factor(gender)), 
                     mod=mod,par.prior=TRUE, prior.plots = F)
#"Removing batch effects and using surrogate variables 
#in differential expression analysis have been shown to 
#reduce dependence, stabilize error rate estimates, and 
#improve reproducibility"
#plotting the "residuals" for one sample
plot(1:nrow(data.matrix),data.matrix[,1]-combat_data[,1],
     xlab="Alz samples",ylab="Corrected difference",col="red")
#plotting the min and max for all samples
plot(uncorrected[1,],col="red",ylim=c(-10,10),main="Uncorrected max and min values for samples")
lines(uncorrected[nrow(uncorrected),],col="green",type="p")
plot(corrected[1,],col="red",ylim=c(-10,10),main="Corrected max and min values for samples")
lines(corrected[nrow(corrected),],col="green",type="p")
#map ids to ensembl genes
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
map<-getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id"),filters="affy_hg_u133_plus_2",values=probesets,ensembl)
#make ABA and Alz matrices compatible for combining
load("~/Google\ Drive/BINF\ 4006\ Project\ /ABAprocessing.RDa")
#uncorrected Alz
dup<-data.frame(data.matrix)
red_dup<-dup[which(probesets %in% map$affy_hg_u133_plus_2),]
genes<-map$ensembl_gene_id[which(probesets %in% map$affy_hg_u133_plus_2)]
red_dup$egenes<-genes
agg<-aggregate(red_dup[,-ncol(red_dup)],by=list(egenes=red_dup$egenes),FUN=mean)
intersected<-intersect(colnames(mat),agg$egenes)
red_agg<-t(agg[which(agg$egenes %in% intersected),])
colnames(red_agg)<-intersected
alzuncorrected<-apply(red_agg[-c(1),],2,as.numeric)
#corrected Alz
dup<-data.frame(fsv1$new) #corrected
red_dup<-dup[which(probesets %in% map$affy_hg_u133_plus_2),]
genes<-map$ensembl_gene_id[which(probesets %in% map$affy_hg_u133_plus_2)]
red_dup$egenes<-genes
agg<-aggregate(red_dup[,-ncol(red_dup)],by=list(egenes=red_dup$egenes),FUN=mean)
intersected<-intersect(colnames(mat),agg$egenes)
red_agg<-t(agg[which(agg$egenes %in% intersected),])
colnames(red_agg)<-intersected
alzfsv1<-apply(red_agg[-c(1),],2,as.numeric)
dup<-data.frame(fsv2$new) #corrected
red_dup<-dup[which(probesets %in% map$affy_hg_u133_plus_2),]
genes<-map$ensembl_gene_id[which(probesets %in% map$affy_hg_u133_plus_2)]
red_dup$egenes<-genes
agg<-aggregate(red_dup[,-ncol(red_dup)],by=list(egenes=red_dup$egenes),FUN=mean)
intersected<-intersect(colnames(mat),agg$egenes)
red_agg<-t(agg[which(agg$egenes %in% intersected),])
colnames(red_agg)<-intersected
alzfsv2<-apply(red_agg[-c(1),],2,as.numeric)
dup<-data.frame(combat_data) #corrected
red_dup<-dup[which(probesets %in% map$affy_hg_u133_plus_2),]
genes<-map$ensembl_gene_id[which(probesets %in% map$affy_hg_u133_plus_2)]
red_dup$egenes<-genes
agg<-aggregate(red_dup[,-ncol(red_dup)],by=list(egenes=red_dup$egenes),FUN=mean)
intersected<-intersect(colnames(mat),agg$egenes)
red_agg<-t(agg[which(agg$egenes %in% intersected),])
colnames(red_agg)<-intersected
alzcorrected<-apply(red_agg[-c(1),],2,as.numeric)
#subsetting intersected genes in ABAData
aba<-apply(mat[,which(colnames(mat) %in% intersected)],2,as.numeric)
#combing ABAData and Alz data
combwcorrected<-rbind(alzcorrected,aba)
combwuncorrected<-rbind(alzuncorrected,aba)
#compute prcomps
pralzuncorrected<-prcomp(alzuncorrected,center = T,scale. = T)
pralzcorrected<-prcomp(alzcorrected,center = T,scale. = T)
pralzfsv1<-prcomp(alzfsv1,center = T,scale. = T)
pralzfsv2<-prcomp(alzfsv2,center = T,scale. = T)
praba<-prcomp(aba,center=T,scale.=T)
prcombwuncorrected<-prcomp(combwuncorrected,center = T,scale. = T)
prcombwcorrected<-prcomp(combwcorrected,center = T,scale. = T)
#plot
tmp<-pralzfsv2
autoplot(tmp,data=alzpheno,shape="proc_regions",
         colour="dz",size=1)+
  #scale_colour_manual(values=cols)+
  ggtitle("Spatial map of aba")+
  theme(
    title=element_text(face=2,size=10),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=8),
    legend.text=element_text(face=2,size=8)
  )
#####save image#####
save.image("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/Alzheimers.RDa")
