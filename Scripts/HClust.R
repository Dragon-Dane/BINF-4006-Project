#PURPOSE: Perform hierarchical clustering on ABA data to see how well expression
#correctly classifies structures based on their embryonic origin

#####SETWD#####
setwd("~/Google\ Drive/BINF\ 4006\ Project\ ")
#####load libraries#####
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")
library(genefilter)
library(ggfortify)
library(ggplot2)
library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#####load data#####
load("~/Documents/Columbia/Courses/TRANSLATIONAL_BIOINFORMATICS/project/RDa/HClust.RDa")
#####source functions#####
source("~/Google\ Drive/BINF\ 4006\ Project\ /Scripts/Functions.R")
#####hierarchical clustering#####
m<-mat
d<-dist(m)
hc<-hclust(d)
plot(hc)
cut<-length(unique(ontology$Structure.Hierarchy.Tier.3))
cutTier3<-cutree(hc,cut)
#####stacked barplot plot of clustering accuracy#####
tab<-table(cutTier3,ontology$Structure.Hierarchy.Tier.3)
s<-apply(tab,2,sum)
frac<-NULL;for(i in 1:length(s)){frac<-cbind(frac,tab[,i]/s[i])}
cols=rainbow(length(unique(ontology$Structure.Hierarchy.Tier.3)))
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Results/clustering_accuracy.pdf")
bp<-barplot(frac,beside=F,col=cols,xlab="Tree Branches",
        ylab="Structure Proportion")
text(bp,rep(1.03,7),s)
text(4,1.1,"number of structures in tree branches")
legend(0,-.03,legend=names(s),xpd=T,text.width=1.8,inset=.1,ncol=4,horiz=F,col=cols,pch=15,cex=.75,bty="n")
dev.off()
#####heatmap of top varying pathways#####
tmp<-t(path_struc_mat)
topVarPaths <- head(order(-rowVars(tmp)),25)
hm<-tmp[topVarPaths,]
hm<-hm-colMeans(hm)
rownames(hm)<-strucCov(rownames(hm))
cols=rainbow(length(unique(ontology$Structure.Hierarchy.Tier.3)))
sidecols <- cols[ ontology$Structure.Hierarchy.Tier.3 ]
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Results/heatmap_topVaryingPathways.pdf",width=15,height=10)
par(cex.main=0.9)
heatmap.2(hm, trace="none", col=bluered(200), 
          ColSideColors=sidecols,labRow=rownames(hm), 
          mar=c(10,2), scale="row",labCol=F,
          main="Top Varying Pathways",
          cexRow=1,margins=c(4,27))
legend(.85,1,legend=names(s),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols,pch=15,cex=.75,bty="n")
dev.off()
#####heatmap of top non-varying pathways#####
tmp<-t(path_struc_mat)
topVarPaths <- head(order(-rowVars(tmp),decreasing=T),25)
hm<-tmp[topVarPaths,]
hm<-hm-colMeans(hm)
rownames(hm)<-strucCov(rownames(hm))
cols=rainbow(length(unique(ontology$Structure.Hierarchy.Tier.3)))
sidecols <- cols[ ontology$Structure.Hierarchy.Tier.3 ]
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Results/heatmap_topNonVaryingPathways.pdf",width=15,height=10)
par(cex.main=0.9)
heatmap.2(hm, trace="none", col=bluered(200),breaks=20, 
          ColSideColors=sidecols,labRow=rownames(hm), 
          mar=c(10,2), scale="row",labCol=F,
          main="Top NonVarying Pathways",
          cexRow=1,margins=c(4,27))
legend(.85,1,legend=names(s),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols,pch=15,cex=.75,bty="n")
dev.off()
#####heatmap of top varying genes#####
tmp<-t(mat)
topVarGenes <- head(order(-rowVars(tmp)),n=500)
hm_tmp<-tmp[topVarGenes,]
hm<-hm_tmp-colMeans(hm_tmp) #subtracts 
#hm_new<-NULL;for(i in 1:414){hm_new<-cbind(hm_new,hm_tmp[,i]-mean(hm_tmp[,i]))};colnames(hm_new)<-colnames(hm_tmp)
ot<-ontology$Structure.Hierarchy.Tier.3
cols=rainbow(length(unique(ot)))
sidecols <- cols[ ot ]
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="complete")}
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Results/heatmap_topVaryingGenes.pdf",width=15,height=10)
par(cex.main=1,font=2)
heatmap.3(hm, trace="none", col=bluered(200),distfun=mydist,
          ColSideColors=matrix(sidecols),labRow=F,dendrogram="column", 
          mar=c(10,2), scale="row",labCol=F,hclustfun=myclust,
          main="",
          cexRow=1)
legend(.9,1.01,legend=levels(ot),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols,pch=15,cex=.75,bty="n")
dev.off()
#####heatmap of top varying genes MET#####
tmp<-t(strucMat(3,"MET_Metencephalon"))
topVarGenes <- head(order(-rowVars(tmp)),n=500)
hm_tmp<-tmp[topVarGenes,]
hm<-hm_tmp-colMeans(hm_tmp) #subtracts 
ot1<-ontology[ontology$Structure.Hierarchy.Tier.3 %in% "MET_Metencephalon","Structure.Hierarchy.Tier.5"]
ot2<-ontology[ontology$Structure.Hierarchy.Tier.3 %in% "MET_Metencephalon","Structure.Hierarchy.Tier.4"]
ot3<-ontology[ontology$Structure.Hierarchy.Tier.3 %in% "MET_Metencephalon","Structure.Hierarchy.Tier.6"]
cols1=rainbow(length(unique(ot1)))
cols2=c("black","grey")
cols3=heat.colors(length(unique(ot3)))
sidecols1 <- cols1[ droplevels(ot1) ]
sidecols2 <- cols2[ droplevels(ot2) ]
sidecols3 <- cols3[ droplevels(ot3) ]
met_clab<-cbind(" "=sidecols1," "=sidecols2)
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="complete")}
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Results/heatmap_topVaryingGenes_met.pdf",width=15,height=10)
par(cex.main=1,font=2)
heatmap.3(hm, trace="none", col=bluered(200),distfun=mydist,
          hclustfun = myclust,ColSideColorsSize=4,
          ColSideColors=met_clab,labRow=F,dendrogram="column", 
          mar=c(7,2), scale="row",labCol=F,cexRow = 1.5,
          main="")
legend(.25,1.02,legend=levels(droplevels(ot2)),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols2,pch=15,cex=1.2,bty="n")
legend(.4,1.02,legend=levels(droplevels(ot1)),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols1,pch=15,cex=1.2,bty="n")
#legend(0,.6,legend=levels(droplevels(ot3)),xpd=T,text.width=1.8,inset=.1,horiz=F,col=cols3,pch=15,cex=1,bty="n")
dev.off()
#####heatmap.3 tutorial#####
#https://www.biostars.org/p/18211/

#Set a working directory for output files
setwd("~/Documents/Columbia/Courses/TRANSLATIONAL_BIOINFORMATICS/project")

#Create a fake dataset for demonstration purposes
prob_matrix=replicate(100, rnorm(20))
drug_names=paste("drug",letters[1:20],sep="_")
patient_ids=paste("patient",c(1:100),sep="_")
rownames(prob_matrix)=drug_names
colnames(prob_matrix)=patient_ids

#Create fake color side bars
drugclass_colors=sample(c("darkorchid","darkred"), length(drug_names), replace = TRUE, prob = NULL)
drugcategory_colors=sample(c("green","darkgreen"), length(drug_names), replace = TRUE, prob = NULL)
subtype_colors=sample(c("red","blue","cyan","pink","yellow","green"), length(patient_ids), replace = TRUE, prob = NULL)
Mcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
Ncolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
Tcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
HER2colors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
PRcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
ERcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
rlab=t(cbind(drugclass_colors,drugcategory_colors))
clab=cbind(subtype_colors,Mcolors,Ncolors,Tcolors,HER2colors,PRcolors,ERcolors)
rownames(rlab)=c("Class","Category")
colnames(clab)=c("Subtype","M","N","T","HER2","PR","ER")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Create heatmap using custom heatmap.3 source code loaded above
pdf(file="heatmap3_example.pdf")
main_title="Drug Response Predictions"
par(cex.main=1)
heatmap.3(prob_matrix, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12),
          Rowv=TRUE, Colv=TRUE, ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
          density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=drug_names, cexRow=1, col=rev(heat.colors(75)),
          ColSideColorsSize=7, RowSideColorsSize=2, KeyValueName="Prob. Response")
legend("topright",legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
       fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

#Example to show that it now also works with just a single column or single row
matrix <- matrix(1:100, byrow=T, nrow=10)
column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
column_annotation <- as.matrix(column_annotation)
colnames(column_annotation) <- c("Variable X")

row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
row_annotation <- as.matrix(t(row_annotation))
rownames(row_annotation) <- c("Variable Y")

heatmap.3(matrix, hclustfun=myclust, distfun=mydist, RowSideColors=row_annotation, ColSideColors=column_annotation)
#####kmeans clustering#####
set.seed(20)
fit <- kmeans(d,7,nstart=20)
table(fit$cluster,ontology$Structure.Hierarchy.Tier.3)
fit$cluster <- as.factor(fit$cluster)
pr<-prcomp(m,center=T,scale=T)
ggplot(pr$x, 
       aes(PC1,PC2,
           colour = fit$cluster)) +
  #scale_shape_manual(values=1:nlevels(ontology$Structure.Hierarchy.Tier.3)) +
  geom_point()
#####save data#####
save.image("~/Documents/Columbia/Courses/TRANSLATIONAL_BIOINFORMATICS/project/RDa/HClust.RDa")
