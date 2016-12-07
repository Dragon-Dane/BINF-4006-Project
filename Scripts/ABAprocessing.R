#PURPOSE: T.B. project-downloading and processing ABA

#####Download and load packages#####
source("https://bioconductor.org/biocLite.R")
biocLite("ABAData")
install.packages("data.table")
require(ABAData)
require(data.table)
#####DataProc#####
data(dataset_adult)
#need to make df of form rows are genes and columns are structure and (i,j) are signal values
DT<-as.data.table(dataset_adult) #Read more at: http://scl.io/ew36cZjC#gs.s_Jyvwo
#making list of ensembl gene-signal values for each structure 
#list of ~16000x1 matrices
#creates duplicates rows because
df<-NULL
system.time(
  df<-lapply(unique(dataset_adult$structure),
             function(x){
               tmp<-DT[structure==x,signal,by=ensembl_gene_id];
               tmp[["signal"]]
             }
  ))
#get gene ids for aggregation of df2
tmp<-DT[structure==unique(dataset_adult$structure)[1],signal,by=ensembl_gene_id];rows<-tmp[["ensembl_gene_id"]]
#now combing the indivdual 'columns' from previous list and combining them together
#combining 16000x1 matrices for 16000x414 matrix where the rows are genes and columns are structures
df2<-NULL
system.time(
  for(i in 1:length(df)){
    df2<-cbind(df2,df[[i]])
  }
)
#really can just exclude duplicates but whatever this works too
system.time(pca.df<-aggregate(df2,by=list(rows),FUN=mean))
rownames(pca.df)<-pca.df$Group.1
mat<-t(pca.df[-c(1)])
#####save image#####
save.image("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/ABAprocessing.RDa")