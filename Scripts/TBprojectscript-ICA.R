#!/usr/bin/Rscript

#####ICA#####
##run on grand central qsub -q res.q -l h_vmem=200g -N mine -j y -o /gpfs/commons/home/giangrecon-420/mine /gpfs/commons/home/giangrecon-420/mine/TBprojectscript-ICA.R
cat("Downloading fastICA...","\n")
#source("https://bioconductor.org/biocLite.R")
install.packages("fastICA",repos='http://cran.us.r-project.org')
library(fastICA)
cat("Loading ABA data...","\n")
load("~/mine/ABAprocessing.RDa")
cat("Computing PCs...\n")
pc.cor<-prcomp(mat,cor=T)
cat("Computing ICs...\n")
ica <- fastICA(mat, n.comp=11)
cat("Saving image...")
save.image("ica.RDa")
cat("Done")
S <- ica$S
dimnames(S) <- list(dimnames(mat)[[1]],paste("Cmp.",1:11,sep=""))
A <- ica$A
plot(ica$X %*% ica$K, main = "PCA components")
plot(S, main = "ICA components")
biplot(S[,1:2],A[,1:2])
