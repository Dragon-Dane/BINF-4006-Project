#PURPOSE: Compute differential expression for parkinsons and schizophrenia datasets.
#Need set of DEGs computed with SVs.

#####SETWD#####
setwd("~/Google\ Drive/BINF\ 4006\ Project\ ")
#####Load libraries#####
library(limma)
library(sva)
library(GEOquery)
library(biomaRt)
library(VennDiagram)
#####load files#####
subgenes<-read.table("Scripts/union_probes_aps.txt",stringsAsFactors = F)[,1]
#####load data#####
load("~/Documents/Columbia/Courses/TRANSLATIONAL_BIOINFORMATICS/project/RDa/DiffExpr.RDa")
#####Schizophrenia-alt#####
sz_alt_gse_set_name <- "GSE21138";
sz_alt_gsem <- getGEO(sz_alt_gse_set_name,GSEMatrix = T);
sz_alt_pheno<-pData(phenoData(sz_alt_gsem[[1]]))
sz_alt_gsenom <- getGEO(sz_alt_gse_set_name,GSEMatrix = F);
names(GPLList(sz_alt_gsenom))
#get platform info
gpl570<-getGEO(names(GPLList(sz_alt_gsenom)))
#potential batch effect contributing variables:
names(which(apply(sz_alt_pheno,2,function(x){length(table(x))>1})))
#making dataset
sz_alt_map<-getBM(attributes=c("illumina_humanref_8_v2","affy_hg_u133_plus_2"),filters="illumina_humanref_8_v2",values=Table(GPLList(sz_alt_gsenom)[[1]])$ID,ensembl)
sz_alt_probesets <- Table(GPLList(sz_alt_gsenom)[[1]])[which(Table(GPLList(sz_alt_gsenom)[[1]])$ID %in% subgenes),"ID"]
sz_alt_tmp1 <- do.call('cbind',lapply(GSMList(sz_alt_gsenom),function(x){
  tab <- Table(x)
  mymatch <- match(sz_alt_probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
sz_alt_tmp2 <- apply(sz_alt_tmp1,2,function(x) {as.numeric(as.character(x))})
sz_alt_tmp3 <- log2(sz_alt_tmp2)
rownames(sz_alt_tmp3)<-sz_alt_probesets
sz_alt_data.matrix<-na.omit(sz_alt_tmp3)
#map ids to ensembl genes
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
sz_alt_map<-getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id"),filters="affy_hg_u133_plus_2",values=sz_alt_probesets,ensembl)

dz<-unlist(lapply(strsplit(as.character(sz_alt_pheno$characteristics_ch1.1)," "),function(x){x[length(x)]}))
dz[grep("illness",dz,fixed=F)]<-"schizophrenia"
sz_alt_mod <- model.matrix(~dz);
sz_alt_mod0 = model.matrix(~1,data=sz_alt_pheno)
#w/o SVs
# Fit linear model with lmFit and eBayes
sz_alt_wosvs_fit <- lmFit(sz_alt_data.matrix,sz_alt_mod);
sz_alt_wosvs_efit <- eBayes(sz_alt_wosvs_fit);
# Get differentially expressed genes with topTable (DEG) for schizophrenia
sz_alt_wosvs_diff_gene_table <- topTable(sz_alt_wosvs_efit, number = nrow(sz_alt_data.matrix), p.value = 1,coef=2);
hist(sz_alt_wosvs_diff_gene_table$P.Value)
sz_alt_wosvs_pi1<-1-pi0est(sz_alt_wosvs_diff_gene_table$P.Value)$pi0
sz_alt_wosvs_diff_exp_genes <- head(rownames(sz_alt_wosvs_diff_gene_table),n=round(sz_alt_wosvs_pi1*length(rownames(sz_alt_wosvs_diff_gene_table)),0));
#w/ SVs
# get SVs if any
sz_alt_n.sv = num.sv(sz_alt_data.matrix,sz_alt_mod,method="leek")
sz_alt_mod0 = model.matrix(~characteristics_ch1.2,data=sz_alt_pheno)
sz_alt_svobj = sva(sz_alt_mat,sz_alt_mod,sz_alt_mod0,n.sv=sz_alt_n.sv) 
#no significant SVs
#####Schizophrenia#####
## Aim of step: download GSE12649 from GEO using functions in GEOquery package. 
# download GDS dataset first(functions you might need: "getGEO")
sz_gse_set_name <- "GSE12649";
sz_gsem <- getGEO(sz_gse_set_name,GSEMatrix = T);
sz_pheno<-pData(phenoData(sz_gsem[[1]]))
sz_gsenom <- getGEO(sz_gse_set_name,GSEMatrix = F);
names(GPLList(sz_gsenom))
#get platform info
gpl96<-getGEO(names(GPLList(sz_gsenom)))
#potential batch effect contributing variables:
names(which(apply(sz_pheno,2,function(x){length(table(x))>1})))
#making dataset
sz_probesets <- Table(GPLList(sz_gsenom)[[1]])[which(Table(GPLList(sz_gsenom)[[1]])$ID %in% subgenes),"ID"]
sz_tmp1 <- do.call('cbind',lapply(GSMList(sz_gsenom),function(x){
  tab <- Table(x)
  mymatch <- match(sz_probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
sz_tmp2 <- apply(sz_tmp1,2,function(x) {as.numeric(as.character(x))})
sz_data.matrix <- log2(sz_tmp2)
rownames(sz_data.matrix)<-sz_probesets
#remove bipolar disorder patients for now
sz_sub_data.matrix<-tmp3[,-which(sz_pheno$source_name_ch1=="bipolar disorder")]
sz_sub_pheno<-sz_pheno[-which(sz_pheno$source_name_ch1=="bipolar disorder"),]
#map ids to ensembl genes
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
sz_map<-getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id"),filters="affy_hg_u133_plus_2",values=sz_probesets,ensembl)
## Aim of step: Obtain the differentially expressed genes between the poorly and highly metastatic samples using functions in limma package
# You might need the following function in this step: "lmFit","eBayes","topTable" (pay attention to the capital letters, substituion with lower letters might output different results). Read the instructions of these functions and learn how to use them. You can complete this step by the following 3 sub-steps. 
# obtain the experimental design of GDS2865 (which sample is tumor, which sample is cancer).
#w/o SVs --> no gender, age etc reported in sz_pheno
sz_mod <- model.matrix(~source_name_ch1,data=sz_pheno);
sz_mod0 = model.matrix(~1,data=sz_pheno)
# Fit linear model with lmFit and eBayes
sz_wosvs_fit <- lmFit(sz_data.matrix,sz_mod);
sz_wosvs_efit <- eBayes(sz_wosvs_fit);
# Get differentially expressed genes with topTable (DEG) for schizophrenia
sz_wosvs_diff_gene_table <- topTable(sz_wosvs_efit, number = nrow(sz_data.matrix), p.value = 1,coef=3);
hist(sz_wosvs_diff_gene_table$P.Value)
sz_wosvs_pi1<-1-pi0est(sz_wosvs_diff_gene_table$P.Value)$pi0
sz_wosvs_diff_exp_genes <- head(rownames(sz_wosvs_diff_gene_table),n=round(sz_wosvs_pi1*length(rownames(sz_wosvs_diff_gene_table)),0));
# Get differentially expressed genes with topTable (DEG) for Bipolar
bp_wosvs_diff_gene_table <- topTable(sz_wosvs_efit, number = nrow(sz_data.matrix), p.value = 1,coef=2);
hist(bp_wosvs_diff_gene_table$P.Value)
bp_wosvs_pi1<-1-pi0est(bp_wosvs_diff_gene_table$P.Value)$pi0
bp_wosvs_diff_exp_genes <- head(rownames(bp_wosvs_diff_gene_table),n=round(bp_wosvs_pi1*length(rownames(bp_wosvs_diff_gene_table)),0));
#####Parkinsons#####
## Aim of step: download GSE28894 from GEO using functions in GEOquery package. 
# download GDS dataset first(functions you might need: "getGEO")
pk_gse_set_name <- "GSE28894";
pk_gsem <- getGEO(pk_gse_set_name,GSEMatrix = T);
pk_pheno<-pData(phenoData(pk_gsem[[1]]))
pk_gsenom <- getGEO(pk_gse_set_name,GSEMatrix = F);
#get platform info
gpl6104<-getGEO(names(GPLList(pk_gsenom)))
#potential batch effect contributing variables:
names(which(apply(pk_pheno,2,function(x){length(table(x))>1})))
#making dataset
pk_probesets <- pk_map[which(pk_map$affy_hg_u133_plus_2 %in% subgenes),1]
pk_tmp1 <- do.call('cbind',lapply(GSMList(pk_gsenom),function(x){
  tab <- Table(x)
  mymatch <- match(pk_probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
pk_tmp2 <- apply(pk_tmp1,2,function(x) {as.numeric(as.character(x))})
pk_tmp3 <- log2(pk_tmp2)
rownames(pk_tmp3)<-pk_probesets
pk_data.matrix <- na.omit(pk_tmp3)
#map ids to ensembl genes
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
pk_map<-getBM(attributes=c("illumina_humanref_8_v3","affy_hg_u133_plus_2","ensembl_gene_id"),filters="illumina_humanref_8_v3",values=rownames(pk_data.matrix),ensembl)
## Aim of step: Obtain the differentially expressed genes between the poorly and highly metastatic samples using functions in limma package
# You might need the following function in this step: "lmFit","eBayes","topTable" (pay attention to the capital letters, substituion with lower letters might output different results). Read the instructions of these functions and learn how to use them. You can complete this step by the following 3 sub-steps. 
# obtain the experimental design of GDS2865 (which sample is tumor, which sample is cancer).
pk_pheno$characteristics_ch1.1<-relevel(pk_pheno$characteristics_ch1.1,"brain region: striatum")
pk_mod <- model.matrix(~characteristics_ch1.1*characteristics_ch1.2,data=pk_pheno);
#pk_mod <- model.matrix(~characteristics_ch1.2,data=pk_pheno);
pk_mod0 = model.matrix(~characteristics_ch1.3,data=pk_pheno)
#w/o SVs
# Fit linear model with lmFit and eBayes
pk_wosvs_fit <- lmFit(pk_data.matrix,pk_mod);
pk_wosvs_efit <- eBayes(pk_wosvs_fit);
# Get differentially expressed genes with topTable (DEG)
pk_wosvs_diff_gene_table <- topTable(pk_wosvs_efit, number = nrow(pk_data.matrix), p.value = 1,coef=6);
#pk_wosvs_diff_gene_table <- topTable(pk_wosvs_efit, number = nrow(pk_data.matrix), p.value = 1);
hist(pk_wosvs_diff_gene_table$P.Value)
pk_wosvs_pi1<-1-pi0est(pk_wosvs_diff_gene_table$P.Value)$pi0
pk_wosvs_diff_exp_genes <- head(pk_wosvs_diff_gene_table$ID,n=round(pk_wosvs_pi1*length(rownames(pk_wosvs_diff_gene_table)),0));
#w/ SVs
batch = pk_pheno$characteristics_ch1.3
pk_modcombat = model.matrix(~1, data=pk_pheno)
pk_combat_edata = ComBat(dat=pk_data.matrix, batch=batch, mod=pk_modcombat, par.prior=TRUE, prior.plots = TRUE) 
pk_mat<-pk_combat_edata
# get SVs if any
pk_n.sv = num.sv(pk_mat,pk_mod,method="leek")
pk_svobj = sva(pk_mat,pk_mod,pk_mod0,n.sv=pk_n.sv) #non-significant SVs
pk_modSv = cbind(pk_mod,pk_svobj$sv)
pk_mod0Sv = cbind(pk_mod0,pk_svobj$sv)
pk_contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,pk_svobj$n.sv-1)),"C2"=c(0,-1,1,rep(0,pk_svobj$n.sv-1)),"C3"=c(-1,0,1,rep(0,pk_svobj$n.sv-1)))
# Fit linear model with lmFit and eBayes
pk_wsvs_fit <- lmFit(pk_mat,pk_modSv);
pk_fitContrasts = contrasts.fit(pk_wsvs_fit,pk_contrast.matrix)
pk_wsvs_efit <- eBayes(pk_wsvs_fit);
# Get differentially expressed genes with topTable (DEG)
pk_wsvs_diff_gene_table <- topTable(pk_wsvs_efit, number = nrow(pk_mat), p.value = 1);
hist(pk_wsvs_diff_gene_table$P.Value)
pk_wsvs_pi1<-1-pi0est(pk_wsvs_diff_gene_table$P.Value)$pi0
pk_wsvs_diff_exp_genes <- head(pk_wsvs_diff_gene_table$ID,n=round(pk_wsvs_pi1*length(rownames(pk_wsvs_diff_gene_table)),0));
#####Alzheimers#####
#get datasets and pheno data
alz_gsem=getGEO("GSE48350",GSEMatrix = T)
alz_pheno<-pData(phenoData(alz_gsem[[1]]))
alz_gsenom= getGEO("GSE48350",GSEMatrix = F);
#get platform info
gpl6104<-getGEO(names(GPLList(alz_gsenom)))
#potential batch effect contributing variables:
names(which(apply(alz_pheno,2,function(x){length(table(x))>1})))
#map ids to ensembl genes
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
alz_probesets <- Table(GPLList(alz_gsenom)[[1]])[which(Table(GPLList(alz_gsenom)[[1]])$ID %in% subgenes),"ID"]
alz_map<-getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id"),filters="affy_hg_u133_plus_2",values=alz_probesets,ensembl)
#making dataset
alz_data.matrix <- do.call('cbind',lapply(GSMList(alz_gsenom),function(x){
  tab <- Table(x)
  mymatch <- match(alz_probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
alz_data.matrix <- apply(alz_data.matrix,2,function(x) {as.numeric(as.character(x))})
alz_data.matrix <- log2(alz_data.matrix)
rownames(alz_data.matrix)<-alz_probesets
#get sample batch variables
agetmp<-unlist(lapply(strsplit(as.character(alz_pheno$title),"_"),function(x){x[3]}))
age<-gsub("yrs","",agetmp)
gender<-ifelse(grepl("_female",alz_pheno$title),"female","male")
regions<-unlist(lapply(strsplit(as.character(alz_pheno$title),"_"),function(x){x[1]}))
proc_regions<-gsub("-","",gsub(" ","",tolower(regions)))
dz<-ifelse(grepl("_AD_",alz_pheno$title),"AD","ctrl")
stage<-unlist(lapply(strsplit(as.character(alz_pheno$characteristics_ch1.3)," "),function(x){x[3]}))
apoe<-unlist(lapply(strsplit(as.character(alz_pheno$characteristics_ch1.4)," "),function(x){x[3]}))
alzpheno<-data.frame(
  gender,
  proc_regions,
  dz,
  age,
  stage,
  apoe
)
#making model matrices for batch correction with sva
alz_mod<-model.matrix(~dz, data=alzpheno)
alz_mod0 = model.matrix(~gender+proc_regions,data=alzpheno) 
#w/o SVs
# Fit linear model with lmFit and eBayes
alz_wosvs_fit <- lmFit(alz_data.matrix,alz_mod);
alz_wosvs_efit <- eBayes(alz_wosvs_fit);
# Get differentially expressed genes with topTable (DEG)
alz_wosvs_diff_gene_table <- topTable(alz_wosvs_efit, number = nrow(alz_data.matrix), p.value = 1);
hist(alz_wosvs_diff_gene_table$P.Value)
alz_wosvs_pi1<-1-pi0est(alz_wosvs_diff_gene_table$P.Value)$pi0
alz_wosvs_diff_exp_genes <- head(rownames(alz_wosvs_diff_gene_table),n=round(alz_wosvs_pi1*length(rownames(alz_wosvs_diff_gene_table)),0));
#w/ SVs
alz_mod <- model.matrix(~dz,data=alzpheno);
alz_mod0 = model.matrix(~1,data=alzpheno)
batch = alzpheno$gender
alz_modcombat = model.matrix(~1, data=alzpheno)
alz_combat_edata = ComBat(dat=alz_data.matrix, batch=batch, mod=alz_modcombat, par.prior=TRUE, prior.plots = TRUE) 
alz_mat<-alz_combat_edata
# get SVs if any
alz_n.sv = num.sv(alz_mat,alz_mod,method="leek")
alz_svobj = sva(alz_mat,alz_mod,alz_mod0,n.sv=alz_n.sv) #non-significant SVs
alz_modSv = cbind(alz_mod,alz_svobj$sv)
alz_mod0Sv = cbind(alz_mod0,alz_svobj$sv)
alz_contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,alz_svobj$n.sv-1)),"C2"=c(0,-1,1,rep(0,alz_svobj$n.sv-1)),"C3"=c(-1,0,1,rep(0,alz_svobj$n.sv-1)))
# Fit linear model with lmFit and eBayes
alz_wsvs_fit <- lmFit(alz_mat,alz_modSv);
alz_fitContrasts = contrasts.fit(alz_wsvs_fit,alz_contrast.matrix)
alz_wsvs_efit <- eBayes(alz_wsvs_fit);
# Get differentially expressed genes with topTable (DEG)
alz_wsvs_diff_gene_table <- topTable(alz_wsvs_efit, number = nrow(alz_mat), p.value = 1);
hist(alz_wsvs_diff_gene_table$P.Value)
alz_wsvs_pi1<-1-pi0est(alz_wsvs_diff_gene_table$P.Value)$pi0
alz_wsvs_diff_exp_genes <- head(rownames(alz_wsvs_diff_gene_table),n=round(alz_wsvs_pi1*length(rownames(alz_wsvs_diff_gene_table)),0));
#####Cerebellar Cortex vs Non Cerebellar Cortex#####
aba_map<-data.frame(colnames(mat))
aba_data.matrix<-t(mat)
aba_pheno<-data.frame(
  "Cerebellum"=as.integer(grepl("Cx_Cerebral Cortex",ontology$Structure.Hierarchy.Tier.4,fixed=F))
)
aba_mod <- model.matrix(~Cerebellum,data=aba_pheno);
#w/o SVs
# Fit linear model with lmFit and eBayes
aba_wosvs_fit <- lmFit(aba_data.matrix,aba_mod);
aba_wosvs_efit <- eBayes(aba_wosvs_fit);
# Get differentially expressed genes with topTable (DEG)
aba_wosvs_diff_gene_table <- topTable(aba_wosvs_efit, number = nrow(aba_data.matrix), p.value = 1);
hist(aba_wosvs_diff_gene_table$P.Value)
aba_wosvs_pi1<-1-pi0est(aba_wosvs_diff_gene_table$P.Value)$pi0
aba_wosvs_diff_exp_genes <- head(rownames(aba_wosvs_diff_gene_table),n=round(aba_wosvs_pi1*length(rownames(aba_wosvs_diff_gene_table)),0));
#####converting ids to ensembl#####
groups<-c("alz","pk","sz","aba","alz_svs","bp","sz_alt")
maps<-list("alz"=alz_map,"pk"=pk_map,"sz"=sz_map,"aba"=aba_map,"alz_svs"=alz_map,"bp"=sz_map,"sz_alt"=sz_alt_map)
dgenes<-list("alz"=alz_wosvs_diff_exp_genes,"pk"=pk_wosvs_diff_exp_genes,"sz"=sz_wosvs_diff_exp_genes,"aba"=aba_wosvs_diff_exp_genes,"alz_svs"=alz_wsvs_diff_exp_genes,"bp"=bp_wosvs_diff_exp_genes,"sz_alt"=sz_alt_wosvs_diff_exp_genes)
sub_dgenes<-NULL
for(i in groups){
  sub_dgenes[[i]]<-unique(maps[[i]][which(maps[[i]][,1] %in% dgenes[[i]]),ncol(maps[[i]])])
}
ensembl_subgenes<-unique(getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id"),filters="affy_hg_u133_plus_2",values=subgenes,ensembl)$ensembl_gene_id)
#####DiffExpr gene pval histograms#####
pdf(paste0(getwd(),"/Results/Proposal2_PkandSz_dgenes_pval_hists.pdf"),height=5,width=10)
par(mfrow=c(1,2))
hist(pk_wosvs_diff_gene_table$P.Value,main="Parkinsons disease",xlab="nominal p-values")
hist(sz_wosvs_diff_gene_table$P.Value,main="Schizophrenia",xlab="nominal p-values")
dev.off()
#####Pathway Enrichment#####
## Aim of step: perform patwhay enrichment analysis
# Read in the file and create pathway dict
pathlist<-read.delim("~/Google\ Drive/BINF\ 4006\ Project\ /Ensembl2Reactome_All_Levels.txt",sep="\t",stringsAsFactors=F,header=F)
reactome<-pathlist[which(pathlist[,ncol(pathlist)]=="Homo sapiens"),]
reactome_dict<-makeDictionary(reactome,3,1)
path_names_in_dict<-reactome[which(reactome[,3] %in% names(reactome_dict)),4]
group_list<-list("alz"=sub_dgenes[["alz"]],"pk"=sub_dgenes[["pk"]],
               "sz"=sub_dgenes[["sz"]],"aba"=aba_wosvs_diff_exp_genes,
               "alz_svs"=sub_dgenes[["alz_svs"]],"bp"=sub_dgenes$bp,
               "sz_alt"=sub_dgenes[["sz_alt"]])
enrichments<-NULL
for(i in 1:length(group_list)){
  candidateGenes<-group_list[[i]]
  backgroundGenes<-setdiff(ensembl_subgenes,group_list[[i]])
  enrchtest<-NULL
  for(j in 1:length(reactome_dict)){
    pathGenes<-reactome_dict[[j]]
    p<-fisher.exact(pathGenes,candidateGenes,backgroundGenes)$p.value
    path_id<-names(reactome_dict)[j]
    path_name<-unique(reactome[which(reactome$V3 %in% path_id),4])
    enrchtest<-rbind(enrchtest,c("cGenes"=length(candidateGenes),
                                 "bGenes"=length(backgroundGenes),
                                 "pathGenes"=length(pathGenes),
                                 "path_id"=path_id,"path_name"=path_name,
                                 "pval"=p))
  }
  padj<-p.adjust(as.numeric(as.character(enrchtest[,ncol(enrchtest)])),method="BH")
  enrichments[[i]]<-as.data.frame(enrchtest)
  enrichments[[i]]$padj<-padj
}
names(enrichments)<-names(group_list)
#plot pvals
pdf(paste0(getwd(),"/Results/diff_sets_pval_hists.pdf"),height=7,width=7)
par(mfrow=c(2,2))
for(i in 1:length(enrichments)){
  hist(as.numeric(as.character(enrichments[[i]]$pval)),xlab="pvals",main=names(enrichments)[i],breaks=100)
}
dev.off()
topEnrichments<-NULL
for(i in 1:length(enrichments)){
  tmp<-enrichments[[i]]
  #filtered<-df[which(as.numeric(as.character(df$pathGenes))>5&as.numeric(as.character(df$pathGenes))<500),]
  df<-tmp
  padjs<-as.numeric(as.character(df$padj))
  #hist(pvals,main="",xlab="")
  topEnrichments[[i]]<-head(df[order(padjs),],n=nrow(df[which(df$padj<0.05),]))
}
names(topEnrichments)<-names(group_list)
#####Venn Diagrams#####
s1<-sub_dgenes$pk
s2<-sub_dgenes$sz
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSz_DEGenes.pdf"),width=4,height=4)
draw.pairwise.venn(length(s1),length(s2),length(intersect(s1,s2)),
                   category=c("Parkinsons","Schizophrenia"),col=c("red","blue"),
                   cat.pos = c(190,170),scaled=T)
dev.off()
s1<-topEnrichments$pk$path_id
s2<-topEnrichments$sz$path_id
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSz_EPathIds.pdf"),width=4,height=4)
draw.pairwise.venn(length(s1),length(s2),length(intersect(s1,s2)),
                   category=c("Parkinsons","Schizophrenia"),col=c("red","blue"),
                   cat.pos = c(190,160),scaled=T)
dev.off()
s1<-sub_dgenes$pk
s2<-sub_dgenes$sz
s3<-sub_dgenes$aba
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSzandABA_DEGenes.pdf"),width=5,height=5)
draw.triple.venn(length(s1),length(s2),length(s3),length(intersect(s1,s2)),
                 length(intersect(s2,s3)),length(intersect(s1,s3)),length(intersect(s1,intersect(s2,s3))),
                 category=c("Parkinsons","Schizophrenia","Cerebellar Cortex"),col=c("red","blue","green"),
                 cat.pos = c(330,30,180),scaled=T)
dev.off()
s1<-topEnrichments$pk$path_id
s2<-topEnrichments$sz$path_id
s3<-topEnrichments$aba$path_id
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSzandABA_EPathIds.pdf"),width=5,height=5)
draw.triple.venn(length(s1),length(s2),length(s3),length(intersect(s1,s2)),
                 length(intersect(s2,s3)),length(intersect(s1,s3)),length(intersect(s1,intersect(s2,s3))),
                   category=c("Parkinsons","Schizophrenia","Cerebellar Cortex"),col=c("red","blue","green"),
                   cat.pos = c(330,30,180),scaled=T)
dev.off()

s1<-sub_dgenes$pk
s2<-sub_dgenes$sz
s3<-sub_dgenes$bp
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSzandBP_DEGenes.pdf"),width=5,height=5)
draw.triple.venn(length(s1),length(s2),length(s3),length(intersect(s1,s2)),
                 length(intersect(s2,s3)),length(intersect(s1,s3)),length(intersect(s1,intersect(s2,s3))),
                 category=c("Parkinsons","Schizophrenia","Bipolar Disorder"),col=c("red","blue","green"),
                 cat.pos = c(330,30,180),scaled=T)
dev.off()
s1<-topEnrichments$pk$path_id
s2<-topEnrichments$sz$path_id
s3<-topEnrichments$bp$path_id
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSzandBP_EPathIds.pdf"),width=5,height=5)
draw.triple.venn(length(s1),length(s2),length(s3),length(intersect(s1,s2)),
                 length(intersect(s2,s3)),length(intersect(s1,s3)),length(intersect(s1,intersect(s2,s3))),
                 category=c("Parkinsons","Schizophrenia","Bipolar Disorder"),col=c("red","blue","green"),
                 cat.pos = c(330,30,360),scaled=T)
dev.off()
s1<-sub_dgenes$pk
s2<-sub_dgenes$sz
s3<-sub_dgenes$bp
s4<-sub_dgenes$aba
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSzandBPandABA_DEGenes.pdf"),width=7,height=5)
draw.quad.venn(length(s1),length(s2),length(s3),length(s4),length(intersect(s1,s2)),
                 length(intersect(s1,s3)),length(intersect(s1,s4)),
               length(intersect(s2,s3)),length(intersect(s2,s4)),length(intersect(s3,s4)),
               length(intersect(s1,intersect(s2,s3))),length(intersect(s1,intersect(s2,s4))),
               length(intersect(s1,intersect(s3,s4))),length(intersect(s2,intersect(s3,s4))),
               length(intersect(s1,intersect(s2,intersect(s3,s4)))),
                 category=c("Parkinsons","Schizophrenia","Bipolar Disorder","Cerebellar Cortex"),col=c("red","blue","purple","green"),
                 cat.pos = c(350,7,360,360),scaled=T)
dev.off()
s1<-topEnrichments$pk$path_id
s2<-topEnrichments$sz$path_id
s3<-topEnrichments$bp$path_id
s4<-topEnrichments$aba$path_id
pdf(paste0(getwd(),"/Results/Proposal2_VD_PkandSzandBPandABA_EPathIds.pdf"),width=7,height=5)
draw.quad.venn(length(s1),length(s2),length(s3),length(s4),length(intersect(s1,s2)),
               length(intersect(s1,s3)),length(intersect(s1,s4)),
               length(intersect(s2,s3)),length(intersect(s2,s4)),length(intersect(s3,s4)),
               length(intersect(s1,intersect(s2,s3))),length(intersect(s1,intersect(s2,s4))),
               length(intersect(s1,intersect(s3,s4))),length(intersect(s2,intersect(s3,s4))),
               length(intersect(s1,intersect(s2,intersect(s3,s4)))),
               category=c("Parkinsons","Schizophrenia","Bipolar Disorder","Cerebellar Cortex"),col=c("red","blue","purple","green"),
               cat.pos = c(350,7,360,360),scaled=T)
dev.off()
#####Overlap significance#####
#https://www.r-bloggers.com/hypothesis-testing-on-normally-distributed-data-in-r/
n.test.oneTails <- function(data, mu0, alpha){
  z<-(mean(data)-mu0)/sd(data)
  n.critical <- qnorm(1-alpha/2) #Es alpha 0.05 -> -1.9599 (df=Inf)
  if(z > -n.critical){
    p.value<-1-pnorm(z)
  }
  else{
    p.value<-pnorm(z)
  }
  if(abs(z) >= n.critical){
    h<-"Reject"
  }
  else{
    h<-"Accept"
  }
  list("H0"=h,"n_critical_values"=c(-n.critical,n.critical),
       "z_score"=z,"p_value"=p.value)
}
s1<-sub_dgenes$pk;lens1<-length(s1)
s2<-sub_dgenes$sz;lens2<-length(s2)
obs<-length(intersect(s1,s2))
null<-c()
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  null<-c(null,length(intersect(t1,t2)))
}
hist(null,xlim=c(min(null)-350,max(null)+350));abline(v=obs,col="red",lwd=2)
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
s1<-topEnrichments$pk$path_id;lens1<-length(s1)
s2<-topEnrichments$sz$path_id;lens2<-length(s2)
obs<-length(intersect(s1,s2))
null<-c() #null distribution is half normal...odd
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  null<-c(null,length(intersect(t1,t2)))
}
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
hist(null,xlim=c(min(null)-350,max(null)+350))
abline(v=obs,col="red",lwd=2)
s1<-sub_dgenes$aba;lens1<-length(s1)
s2<-sub_dgenes$sz;lens2<-length(s2)
obs<-length(intersect(s1,s2))
null<-c()
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  null<-c(null,length(intersect(t1,t2)))
}
hist(null,xlim=c(min(null)-350,max(null)+350));abline(v=obs,col="red",lwd=2)
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
s1<-sub_dgenes$pk;lens1<-length(s1)
s2<-sub_dgenes$aba;lens2<-length(s2)
obs<-length(intersect(s1,s2))
null<-c()
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  null<-c(null,length(intersect(t1,t2)))
}
hist(null,xlim=c(min(null)-350,max(null)+350));abline(v=obs,col="red",lwd=2)
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
s1<-topEnrichments$pk$path_id;lens1<-length(s1)
s2<-topEnrichments$aba$path_id;lens2<-length(s2)
obs<-length(intersect(s1,s2))
null<-c() #null distribution is half normal...odd
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  null<-c(null,length(intersect(t1,t2)))
}
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
hist(null,xlim=c(min(null)-350,max(null)+350))
abline(v=obs,col="red",lwd=2)
s1<-topEnrichments$aba$path_id;lens1<-length(s1)
s2<-topEnrichments$sz$path_id;lens2<-length(s2)
obs<-length(intersect(s1,s2))
null<-c() #null distribution is half normal...odd
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  null<-c(null,length(intersect(t1,t2)))
}
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
hist(null,xlim=c(min(null)-350,max(null)+350))
abline(v=obs,col="red",lwd=2)

s1<-sub_dgenes$pk;lens1<-length(s1)
s2<-sub_dgenes$sz;lens2<-length(s2)
s3<-sub_dgenes$aba;lens3<-length(s3)
obs<-length(intersect(s1,intersect(s2,s3)));null<-c()
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  t3<-sample(ensembl_subgenes,lens3)
  null<-c(null,length(intersect(t1,intersect(t2,t3))))
}
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
hist(null,xlim=c(min(null)-350,max(null)+350))
abline(v=obs,col="red",lwd=2)

s1<-topEnrichments$pk$path_id;lens1<-length(s1)
s2<-topEnrichments$sz$path_id;lens2<-length(s2)
s3<-topEnrichments$aba$path_id;lens3<-length(s3)
obs<-length(intersect(s1,intersect(s2,s3)));null<-c()
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  t3<-sample(ensembl_subgenes,lens3)
  null<-c(null,length(intersect(t1,intersect(t2,t3))))
}
options(scipen = -999999)
n.test.oneTails(null,obs,1e-5)
hist(null,xlim=c(min(null)-350,max(null)+350))
abline(v=obs,col="red",lwd=2)

#simulate pvals
s1<-sub_dgenes$pk
lens1<-length(s1)
s2<-sub_dgenes$sz
lens2<-length(s2)
obs<-length(intersect(s1,s2))
null<-c()
for(i in 1:1000){
  t1<-sample(ensembl_subgenes,lens1)
  t2<-sample(ensembl_subgenes,lens2)
  null<-c(null,length(intersect(t1,t2)))
}
p<-c()
z<-c()
for(i in (obs-100):(obs+100)){
  p<-c(p,n.test.oneTails(null,i,0.05)$p_value)
  z<-c(z,n.test.oneTails(null,i,0.05)$z_score)
}
plot(p,pch=19)
#####Top common enriched pathway in PK and SZ#####
s1<-topEnrichments$pk
s2<-topEnrichments$sz
top1<-as.character(s1[which.min(s1[which(intersect(s1$path_id,s2$path_id) %in% s1$path_id),"padj"]),"path_id"])
top2<-as.character(s2[which.min(s2[which(intersect(s1$path_id,s2$path_id) %in% s2$path_id),"padj"]),"path_id"])
#####Top diff genes in PK and SZ and BP-output for PheGenI#####
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
genes<-head(sub_dgenes$pk,n=100)
map<-getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters="ensembl_gene_id",values=genes,ensembl)
cat(map[,2])
#####PheGenI graphs#####
pk_phegeni<-read.table(paste0(getwd(),"/Results/PK_top100_DEGs_PheGenI_GaP.tab"),stringsAsFactors=F,fill=T,sep="\t",header=T)
sz_phegeni<-read.table(paste0(getwd(),"/Results/SZ_top100_DEGs_PheGenI_GaP.tab"),stringsAsFactors=F,fill=T,sep="\t",header=T)
list<-pk_phegeni
pk_tmp_tab<-unlist(lapply(strsplit(list$Disease,";"),paste))
pk_tab<-table(pk_tmp_tab)[order(table(pk_tmp_tab),decreasing = T)]
write.table(data.frame(pk_tab),paste0(getwd(),"/Results/pk_phegeni_table.txt"),quote=F,row.names=F,sep="\t",col.names=F)
list<-sz_phegeni
sz_tmp_tab<-unlist(lapply(strsplit(list$Disease,";"),paste))
sz_tab<-table(sz_tmp_tab)[order(table(sz_tmp_tab),decreasing = T)]
write.table(data.frame(sz_tab),paste0(getwd(),"/Results/sz_phegeni_table.txt"),quote=F,row.names=F,sep="\t",col.names=F)
s1<-gsub(" ","",names(head(sz_tab,n=5)))
s2<-gsub(" ","",names(head(pk_tab,n=5)))
length(intersect(s1,s2))
#####Save image#####
save.image("~/Documents/Columbia/Courses/TRANSLATIONAL_BIOINFORMATICS/project/RDa/DiffExpr.RDa")
