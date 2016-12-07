#PURPOSE: T.B. project-compute PCA for data

#####load Data#####
load("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/ABAprocessing.RDa")
#####basic ABA PCA#####
pr_mat<-prcomp(mat,
           center = T,
           scale. = T)
plot(pr_mat, type = "l")
summary(pr_mat)$importance[,1:2]
autoplot(pr_mat,type="obs")
#####Descriptive statistics####
hist(apply(mat,1,mean),breaks=20,xlab="mean in 414 brain structures")
hist(apply(mat,1,var),breaks=20,xlab="variance in 414 brain structures")
hist(apply(mat,1,sd),breaks=20,xlab="sd in 414 brain structures")
plot(pr_mat, type = "l")
pr_mat.var<-pr$sdev^2
#proportion of variance explained
pve<-pr_mat.var/sum(pr_mat.var)
#plot of above
plot(cumsum(pve ), xlab="Principal Component ", ylab=" Cumulative Proportion of Variance Explained ", ylim=c(0,1) , type="b")
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
#####Visualize and print-to-pdf#####
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Proposal1_PCA_ABAData_414structures_ColoredBySuperstructure3_ShapeBySuperstructure2.pdf",width=10,height=7)
cols<-c("red","blue","violet","green4","darkgray","orange","black")
autoplot(pr_mat,data=ontology,shape="Structure.Hierarchy.Tier.2",
         colour="Structure.Hierarchy.Tier.3",size=3)+
  scale_colour_manual(values=cols)+
  ggtitle("Spatial map of brain structure signatures")+
  theme(
    title=element_text(face=2,size=15),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=15),
    legend.text=element_text(face=2,size=10)
  )
dev.off()
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Proposal1_PCA_ABAData_414structures_ColoredBySuperstructure4_ShapeBySuperstructure3.pdf",width=17,height=12)
cols<-c("red","blue","violet","green4","darkgray","darkgreen","black","mediumpurple","tomato3","darkgoldenrod","salmon4","orchid","skyblue","peachpuff2","darkkhaki","maroon","springgreen4","firebrick","darkorange","turquoise","snow4","magenta2","yellow2","seagreen3","purple","navy")
autoplot(pr_mat,data=ontology,shape="Structure.Hierarchy.Tier.3",
         colour="Structure.Hierarchy.Tier.4",size=4)+
  scale_colour_manual(values=cols)+
  scale_shape_manual(values=1:nlevels(ontology$Structure.Hierarchy.Tier.3)) +
  ggtitle("Distinguishing finer brain structures within major superstructures")+
  theme(
    title=element_text(face=2,size=18),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=20),
    legend.text=element_text(face=2,size=12)
  )
dev.off()
#####subsetting PCA by Metencephalon and cerebellum#####
met_indices<-which(ontology$Structure.Hierarchy.Tier.3 %in% "MET_Metencephalon")
newontology<-ontology[which(ontology$Structure.Hierarchy.Tier.3 %in% "MET_Metencephalon"),]
pca_met<-mat[met_indices,]
pr_met<-prcomp(pca_met,
               center = T,
               scale. = T)
plot(pr_met, type = "l")
summary(pr_met)$importance[,1:2]
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Proposal2_PCA_ABAData_Metencephalon_ColoredBySuperstructure6_ShapeBySuperstructure5.pdf",width=10,height=8)
cols<-c("red","blue","violet","orange1","darkgray","darkgreen","black","mediumpurple","darkgoldenrod","grey50","orchid","skyblue","peachpuff2","magenta2","turquoise","yellow2","seagreen3","navy")
autoplot(pr_met,data=newontology,shape="Structure.Hierarchy.Tier.5",
         colour="Structure.Hierarchy.Tier.6",size=3)+
  scale_colour_manual(values=cols)+
  #geom_text(aes(label=newontology$Structure.Hierarchy.Tier.5),
  #          size=3,
  #          vjust=-0.5)+
  ggtitle("Spatial map of metencephalon signatures")+
  theme(
    title=element_text(face=2,size=15),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=15),
    legend.text=element_text(face=2,size=10)
  )
dev.off()
cc_indices<-which(ontology$Structure.Hierarchy.Tier.5 %in% "CbCx_Cerebellar Cortex")
ccontology<-ontology[which(ontology$Structure.Hierarchy.Tier.5 %in% "CbCx_Cerebellar Cortex"),]
pca_cc<-mat[cc_indices,]
pr_cc<-prcomp(pca_cc,
              center = T,
              scale. = T)
summary(pr_cc)$importance[,1:2]
plot(pr_cc, type = "l")
pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Proposal2_PCA_ABAData_CerebellarCortex_ColoredBySuperstructure8_ShapeBySuperstructure7.pdf",width=10,height=8)
cols<-c("red","blue","violet","green4","darkgray","darkgreen","black","mediumpurple","tomato3","darkgoldenrod","salmon4","darkkhaki","maroon","springgreen4","firebrick","darkorange","turquoise","snow4","magenta2","yellow2","seagreen3","purple","navy")
autoplot(pr_cc,data=ccontology,shape="Structure.Hierarchy.Tier.7",
         colour="Structure.Hierarchy.Tier.8",size=3)+
  scale_colour_manual(values=cols)+
  ggtitle("Spatial map of cerebellar cortex signatures")+
  theme(
    title=element_text(face=2,size=15),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=15),
    legend.text=element_text(face=2,size=10)
  )
dev.off()
#mapping for superstructures to individual structures
str3<-unique(ontology$Structure.Hierarchy.Tier.3)[4]
str3
unique(subset(ontology,Structure.Hierarchy.Tier.3==str3)$Structure.Hierarchy.Tier.6)
#####subsetting PCA by Telencephalon#####
tel_indices<-which(ontology$Structure.Hierarchy.Tier.3 %in% "Tel_Telencephalon")
newontology<-ontology[which(ontology$Structure.Hierarchy.Tier.3 %in% "Tel_Telencephalon"),]
pca_tel<-mat[tel_indices,]
pr_tel<-prcomp(pca_tel,
               center = T,
               scale. = T)
plot(pr_tel, type = "l")
summary(pr_met)$importance[,1:2]
#pdf("~/Google\ Drive/BINF\ 4006\ Project\ /Proposal2_PCA_ABAData_Metencephalon_ColoredBySuperstructure6_ShapeBySuperstructure5.pdf",width=10,height=8)
cols<-c("red","blue","violet","orange1","darkgray","darkgreen","black","mediumpurple","darkgoldenrod","grey50","orchid","skyblue","peachpuff2","magenta2","turquoise","yellow2","seagreen3","navy")
autoplot(pr_tel,data=newontology,shape="Structure.Hierarchy.Tier.4",
         colour="Structure.Hierarchy.Tier.5",size=3)+
  #scale_colour_manual(values=cols)+
  #geom_text(aes(label=newontology$Structure.Hierarchy.Tier.5),
  #          size=3,
  #          vjust=-0.5)+
  ggtitle("Spatial map of telencephalon signatures")+
  theme(
    title=element_text(face=2,size=15),
    axis.text=element_text(face=2),
    axis.title=element_text(face=2,size=15),
    legend.text=element_text(face=2,size=10)
  )
#dev.off()

#the Pons and Cerebellum, though in the same superstructure, have very disparate expression
#####This is wicked cool!#####
#http://planspace.org/2013/02/03/pca-3d-visualization-and-clustering-in-r/
pr_mat<-prcomp(mat,
           center = T,
           scale. = T)
plot3d(pr$x[,1:3])
#####save image#####
save.image("~/Google\ Drive/BINF\ 4006\ Project\ /RDa/PCA.RDa")
