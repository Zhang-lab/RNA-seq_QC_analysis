library("DESeq2")
# load count matrix######################################################################################################################
setwd("/Users/chengl/Desktop/")

Bartolomei=read.table("Bartolomei.txt",header=T,sep="\t")
Dolinoy=read.table("Dolinoy.txt",header=T,sep="\t")
Mutlu=read.table("Mutlu.txt",header=T,sep="\t")

rownames(Bartolomei)=Bartolomei[,1]
Bartolomei=Bartolomei[,-1]
Bartolomei=Bartolomei[,-1]

rownames(Dolinoy)=Dolinoy[,1]
Dolinoy=Dolinoy[,-1]
Dolinoy=Dolinoy[,-1]
Dolinoy=Dolinoy[,-which(colnames(Dolinoy)=="T105c_Lead_F_Liver_5mo.R1")]

rownames(Mutlu)=Mutlu[,1]
Mutlu=Mutlu[,-1]
Mutlu=Mutlu[,-1]

countdata=cbind(Bartolomei,Dolinoy,Mutlu)
# load experiment design#####################################################################################################################
tt=read.table("BartolomeiLab_exp_design.txt",header=T,sep="\t")
Tissue=factor(rep(1,17),label="Liver")
tt=cbind(tt,Tissue)
group1=tt[tt$SAMPLE%in%colnames(Bartolomei),2]
sex1=tt[tt$SAMPLE%in%colnames(Bartolomei),3]
tissue1=tt[tt$SAMPLE%in%colnames(Bartolomei),4]

sex=c()
group=c()
tissue=c()
for(i in 1:dim(Dolinoy)[2]) {
  if(length(grep("M",colnames(Dolinoy)[i],fixed=T))==1) {
    sex=c(sex,"MALE")
  } else {
    sex=c(sex,"FEMALE")
  }
  
  if(length(grep("Ctrl",colnames(Dolinoy)[i],fixed=T))==1) {
    group=c(group,"Ctrl")
  } else if(length(grep("Lead",colnames(Dolinoy)[i],fixed=T))==1) {
    group=c(group,"Lead")
  } else{
    group=c(group,"DEHP")
  }

  if(length(grep("Liver",colnames(Dolinoy)[i],fixed=T))==1) {
    tissue=c(tissue,"Liver")
  } else {
    tissue=c(tissue,"Blood")
  }
}
sex2=sex
group2=group
tissue2=tissue

tt=read.table("MutluLab_exp_design.txt",header=T,sep="\t")
Tissue=factor(c(rep(0,24),rep(1,24),rep(2,24)),label=c("Lung","Liver","Heart"))
tt=cbind(tt,Tissue)
group3=tt[tt$Samples%in%colnames(Mutlu),3]
sex3=tt[tt$Samples%in%colnames(Mutlu),4]
tissue3=tt[tt$Samples%in%colnames(Mutlu),5]

group=as.factor(c(as.character(group1),as.character(group2),as.character(group3)))
sex=as.factor(c(as.character(sex1),as.character(sex2),as.character(sex3)))
lab=factor(c(rep(0,ncol(Bartolomei)),rep(1,ncol(Dolinoy)),rep(2,ncol(Mutlu))),label=c("Bartolomei","Dolinoy","Mutlu"))
tissue=as.factor(c(as.character(tissue1),as.character(tissue2),as.character(tissue3)))

colData=data.frame(lab,sex,group,tissue)
rownames(colData)=colnames(countdata)
# create DESeq object and pre-filter#####################################################################################################################
dds=DESeqDataSetFromMatrix(countData=countdata,colData=colData,design=~lab+sex+tissue)
dds=dds[rowSums(fpm(dds,robust=F)>10)>10,]
# transformation#####################################################################################################################
rld=vst(dds,blind=FALSE)
# distance analysis#####################################################################################################################
library("pheatmap")
library("RColorBrewer")
library(grid)

sampleDists=dist(t(assay(rld)))
sampleDistMatrix=as.matrix(sampleDists)
rownames(sampleDistMatrix)=tissue
colnames(sampleDistMatrix)=lab
colors=colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

png("distance_ALL.png",height=3700,width=3700,res=300)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors,main="Heatmap of similarity between ALL samples based on Euclidean distance")
setHook("grid.newpage", NULL, "replace")
grid.text("Index of samples", y=-0.02, gp=gpar(fontsize=12))
grid.text("Index of tissue", x=-0.03, rot=90, gp=gpar(fontsize=12))
dev.off()
# correlation analysis#####################################################################################################################
df=as.data.frame(colData(dds)[,c("sex","lab","tissue")])
png("correlation_ALL.png",height=3700,width=3700,res=300)
pheatmap(cor(assay(rld)),annotation_col=df,show_colnames=F,main="Heatmap of correlation between ALL samples")
dev.off()
# PCA#####################################################################################################################
library(ggplot2)

pcaData=plotPCA(rld,intgroup=c("group","sex","lab","tissue"), returnData = TRUE)
percentVar=round(100*attr(pcaData,"percentVar"))

ggplot(pcaData,aes(x=PC1,y=PC2,size=lab,shape=sex,color=tissue))+
  geom_point()+
  scale_shape_manual(values=c(1,2,16))+
  scale_size_manual(values=c(6,4,8))+
  ggtitle("Principal component analysis with covariate of ALL samples")+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
        text=element_text(size=12,family="Tahoma"),
        axis.title=element_text(face="bold"),
        axis.text.x=element_text(size=10,face="bold"),
        axis.text.y=element_text(size=10,face="bold"),
        legend.text=element_text(size=10,face="bold"),
        legend.title=element_text(size=10,face="bold"))