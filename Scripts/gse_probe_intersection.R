source("https://bioconductor.org/biocLite.R")
library(GEOquery)
library(limma)
library(data.table)
library("biomaRt")
biocLite("annotate")
library("annotate")
biocLite("biomaRt")
# helpful: https://www.biostars.org/p/76097/
esetGSE48350 <- getGEO("GSE48350")
esetGSE48350 <- esetGSE48350[[1]]
expGSE48350 <- exprs(esetGSE48350)
probGSE48350 <- rownames(expGSE48350)
esetGSE28894 <- getGEO("GSE28894")
esetGSE28894 <- esetGSE28894[[1]]
expGSE28894 <- exprs(esetGSE28894)
probGSE28894 <- rownames(expGSE28894)

biocLite("hugene10stprobeset.db")
esetGSE35974 <- getGEO("GSE35974")
esetGSE35974 <- esetGSE35974[[1]]
expGSE35974 <- exprs(esetGSE35974)
prob35974 <- rownames(expGSE35974)
library("hugene10stprobeset.db")
# OUT <- select(hugene10stprobeset.db, prob35974, c("SYMBOL","ENTREZID","GENENAME"))
biocLite("hgu133plus2.db")
biocLite("illuminaHumanv2.db")
library("illuminaHumanv2.db")
OUT_28894 <- select(illuminaHumanv2.db, probGSE28894, c("SYMBOL","ENTREZID","GENENAME"))
dim(OUT_28894)
OUT_28894[1:5,]
# the following still giving issues. Maybe try different .db?
#OUT_35974 <- select(hugene10stprobeset.db, prob35974, c("SYMBOL","ENTREZID","GENENAME"))
#OUT_35974[1:5,]
esetGSE12649 <- getGEO("GSE12649")
esetGSE12649 <- esetGSE12649[[1]]
expGSE12649 <- exprs(esetGSE12649)
probGSE12649 <- rownames(expGSE12649)
library("hgu133a.db")
OUT_48350 <- select(hgu133plus2.db, probGSE48350, c("SYMBOL","ENTREZID","GENENAME"))
OUT_48350[1:5,]
#output column with ENTREZID or column 3 to vector for each gse
entrezid_48350 <- as.numeric(as.vector(OUT_48350[,3]))
entrezid_48350[1:5]
entrezid_28894 <- as.numeric(as.vector(OUT_28894[,3]))
entrezid_28894[1:5]
#create list with the intersection of the entrez IDs of both gse's
#still need to add the remaining gse for schizophrenia
# intersect(entrezid_28894, entrezid_48350)
#inter_28894_48350 <- intersect(entrezid_28894, entrezid_48350)
#length(inter_28894_48350)
#length(entrezid_28894)
#length(entrezid_48350)
# find intersection of all 3 datasets for Alz, Sch, and Park
entrezid_12649 <- as.numeric(as.vector(OUT_12649[,3]))
int_aps <- intersect(intersect(entrezid_12649,entrezid_28894),entrezid_48350)

#now print all the rows of the OUT file from in GSE12649 that have elements from intersection list in the entrez ID column
OUT_12649_int <-OUT_12649[is.element(OUT_12649$ENTREZID,int_aps),]
#this is a list of all the unique entrezids of gse12649 that match the intersection list
#uni_OUT_12649_int <- unique(OUT_12649_int$ENTREZID)
#gives length 12525 which is the same as length of intersection list int_aps
#print out all unique probe IDs in GSE12649 corresponding to ENTREZIDs in the AlzParkSchi intersection
uniprob_OUT_12649_int <- unique(OUT_12649_int$PROBEID)
#repeat for other gse's
OUT_28894_int <-OUT_28894[is.element(OUT_28894$ENTREZID,int_aps),]
uniprob_OUT_28894_int <- unique(OUT_28894_int$PROBEID)
OUT_48350_int <-OUT_48350[is.element(OUT_48350$ENTREZID,int_aps),]
uniprob_OUT_48350_int <- unique(OUT_48350_int$PROBEID)

#create a list that is the union of the three lists of unique probes
union_probes_aps <- union(uniprob_OUT_48350_int,union(uniprob_OUT_12649_int,uniprob_OUT_28894_int))
length(union_probes_aps)
#gives length of 55774 unique probeIDs