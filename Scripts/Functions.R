#PURPOSE: To store functions related to BINF Project.

##########data to be loaded for functions
load("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/Ontology.RDa")
load("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/ABAProcessing.RDa")
load("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/Enrichment.RDa")
####################
makeDictionary<-function(file=file,keys,values){
#'file' is a data frame (generated from loading a file into R typically)
#'keys' is the dictionary keys (unique identifier) to associate to values\
#'values' is an identifier which associates to a unique key
  cat("Making dictionary\n")
  cat("Some keys",head(unique(file[,keys])),"\n")
  cat("Some values",head(unique(file[,values])),"\n")
  dict<-vector(mode="list",length=length(unique(file[,keys])))
  names(dict)<-unique(file[,keys])
  for(i in 1:length(dict)){dict[[i]]<-file[which(file[,keys] %in% names(dict)[i]),values]}
  return(dict)
}
####################
candidateABAGenes <- function(structureRowID){
#extract candidate genes from ABA structures (1 out of 414)
#can be any integer from 1 to 414
  #define structure x gene df
  df<-mat
  #define subsetted df
  arr<-df[structureRowID,]
  names(arr)<-colnames(df)
  cutoffs<-quantile(arr,seq(0,1,.01))[c("1%","99%")]
  candidate_genes<-vector(mode="list",length=length(cutoffs))
  names(candidate_genes)<-names(cutoffs)
  tmp1<-as.numeric(gsub("%","",names(cutoffs)));
  for(i in 1:length(tmp1)){
     if(tmp1[i]<50){
        candidate_genes[[i]]<-arr[which(arr<cutoffs[i])]
     }else{
        candidate_genes[[i]]<-arr[which(arr>cutoffs[i])]
     }
  }
  candidate_genes
}
####################
backgroundABAGenes<-function(cgenes){
#getting background genes mutually
#exclusive from ABA candidate genes
  df<-mat
  colnames(df)[which(!(colnames(df) %in% names(cgenes)))]
}
####################
backgroundDEGenes<-function(DEGs,dataset){
#getting background genes mutually
#exclusive from DE genes
  bgGenes<-dataset
  setdiff(bgGenes,DEGs)
}
####################
#First developing enrichment function
pathway.enrichment.analysis <- function(geneSetList,candidateGenes,backgroundGenes){

}
####################
fisher.exact<-function(pgenes,cgenes,bgenes){
#basic hypergeometric test for
#genes in pathways
  cGenes_in_path<-intersect(candidateGenes,pathGenes)
  cGenes_notin_path<-setdiff(candidateGenes,pathGenes)
  bGenes_in_path<-intersect(backgroundGenes,pathGenes)
  bGenes_notin_path<-setdiff(backgroundGenes,pathGenes)
  cont_table<-matrix(c(length(cGenes_in_path),
                       length(bGenes_in_path),
                       length(cGenes_notin_path),
                       length(bGenes_notin_path)),nrow=2,ncol=2,
                     dimnames=list(c("cGenes","bGenes"),c("in_pathway","not_in_pathway"))
  )
  fisher.test(cont_table)
}
####################
strucCov<-function(x){
#convert reactome path id to name 
#or name to id
  if(grepl("-",x)){
    unique(reactome[which(reactome[,3] %in% x),4])
  }else{
    unique(reactome[which(reactome[,4] %in% x),3])
  }
}
####################
varyingPathwayPCA<-function(pathwayOfInterestId,sgMat,spMat){
#this function is to allow for fast generation of PCA of 
#a given set of structures (all, metencephelon, cerebellum etc. 
#and visualizing a given pathway variation
#'pathwayOfInterestId' has to be pathway id from reactome[,3]
#'sgMat' is the genes by structures matrix of interest
#'spMat' is the pathways by structures matrix of interest
#load Enrichment.RDa first
  
  #assign pathway to variable
  pathwayOfInterestName<-unique(reactome[which(reactome[,3] %in% pathwayOfInterestId),4])
  
  #computing components
  pr<-prcomp(sgMat,
             center = T,
             scale. = T)
  
  #changing colnames (pathway ids) to something interpretable by ggplot
  tmp<-spMat
  strs<-strsplit(colnames(spMat),"-")
  tmp_ids<-unlist(lapply(strs,function(x){paste0(x[1],x[3])}))
  rownames(tmp)<-NULL
  colnames(tmp)<-tmp_ids
  data<-as.data.frame(tmp)
  
  #coloring param definition
  colorBy<-colnames(data)[which(colnames(spMat) %in% pathwayOfInterestId)]
  
  #plotting
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
  
}
####################
strucMat<-function(tier,key){
#function for picking structures for generating
#gene by structure matrix-especially for doing 
#unsupervised learning e.g.
#'tier' is the number of hierarchy in ontology
#'key' is the character value in that tier to subset by
  inds<-which(ontology[,paste0("Structure.Hierarchy.Tier.",tier)] %in% key)
  mat[inds,]
}
####################
pathMat<-function(tier,key){
#function for picking structures for generating
#pathway by structure matrix-especially for doing 
#unsupervised learning e.g.
#'tier' is the number of hierarchy in ontology
#'key' is the character value in that tier to subset by
  #get indices of interest
  inds<-which(ontology[,paste0("Structure.Hierarchy.Tier.",tier)] %in% key)

  #subset pathway matrix
  path_struc_mat[inds,]
}
####################

