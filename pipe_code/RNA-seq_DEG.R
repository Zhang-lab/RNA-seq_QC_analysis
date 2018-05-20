#!/usr/bin/env Rscript
args=commandArgs()

if(grepl("-h",args[6])) {
  #############
  # help page #
  #############
  cat("\n####################################################\n")
  cat("# Rscript for performing DEG analysis using DESeq2 #")
  cat("\n####################################################\n")
  cat("\nVersion: 1.1a (Last modified on 05-17-2018)\n")
  cat("\nUsage: Rscript RNA-seq_DEG.R EXPRESSION_TABLE DESIGN_TABLE CPM_CUTOFF PVAL_CUTOFF L2F_CUTOFF PERFORMGO_LOGIC SPECIES GO_PVAL_CUTOFF\n")
  cat("\nNote: the parameters are positional, which means you have to follow the above order of parameters\n")
  cat("\n\n<parameter explanation>\n")
  cat("\n(1) EXPRESSION_TABLE: expression count table directly generated from `batch_collection_rna.sh` in our RNA-seq pipeline\n")
  cat("\n(2) DESIGN_TABLE: table specified the experiment design which contains at least 2 columns")
  cat("\n    (a) the first column refers to the name of samples, which should be identical to those in EXPRESSION_TABLE")
  cat("\n    (b) the second column refers to the main variate that will be used to make a contrast in DESeq")
  cat("\n    (c) when the second column contains more than two levels, user can give another column called \"contrast\" to specify the contrast")
  cat("\n        use <0> for the reference level, <1> for the contrast level, and <2> for other levels")
  cat("\n        if \"contrast\" was not given, only the first two different levels will be used to make a contrast")
  cat("\n    (d) any more columns are optional and will be considered as covariates to be used for adjustment in the linear model")
  cat("\n        if any covariate contains only one level, it will be removed\n")
  cat("\n(3) CPM_CUTOFF: genes with cpm lower than the given cutoff will be removed in the beginning, such as <1>\n")
  cat("\n(4) PVAL_CUTOFF: genes with adjusted p-value smaller than the given cutoff will be considered as significant genes, such as <0.01>\n")
  cat("\n(5) L2F_CUTOFF: genes with absolute log2 fold-change larger than the given cutoff will be considered as significant genes, such as <1>\n")
  cat("\n(6) PERFORMGO_LOGIC: indicator of performing GO analysis or not, such as <TRUE>/<FALSE>\n")
  cat("\n(7) SPECIES: indicator of data species, including human <Hs>, mouse <Mm> and zebrafish <Dr>\n")
  cat("\n(8) GO_PVAL_CUTOFF: GO terms with adjusted p-value smaller than the given cutoff will be considered as enriched categories, such as <0.05>\n")
  cat("\n\n<input file example>\n")
  cat("\n(1) EXPRESSION_TABLE:\n")
  cat("\n    name    length    sample1    sample2    sample3    sample4    sample5    sample6\n")
  cat("    gene1     1111         12         23         15         22         13         34\n")
  cat("    gene2     2222         11         13         33         34         54         15\n")
  cat("\n(2) DESIGN_TABLE:\n")
  cat("\n    name    treatment    contrast    sex\n")
  cat("    sample1        G1           0      m\n")
  cat("    sample2        G1           0      f\n")
  cat("    sample3        G2           2      m\n")
  cat("    sample4        G2           2      f\n")
  cat("    sample5        G3           1      m\n")
  cat("    sample6        G3           1      f\n")
  cat("\nNote: in above example, we will make contrast `G1 vs G3` and adjust for sex. If \"contrast\" was not used, we will make contrast `G1 vs G2`\n")
  cat("\n#######################\n")
  cat("# Thanks for using!!! #")
  cat("\n#######################\n")
  cat("\n")
} else {
  ###########################
  # load required libraries #
  ###########################
  cat("\nLoading libraries ...\n")
  suppressMessages(library(DESeq2))
  suppressMessages(library(ggplot2))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(pheatmap))
  suppressMessages(library(dplyr))
  suppressMessages(library(GOstats))
  suppressWarnings(suppressMessages(library(biomaRt)))
  suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
  suppressWarnings(suppressMessages(library(org.Mm.eg.db)))
  suppressWarnings(suppressMessages(library(org.Dr.eg.db)))
  suppressWarnings(suppressMessages(library(KEGG.db)))
  suppressWarnings(suppressMessages(library(gage)))
  suppressMessages(library(DBI))
  
  ##########################################
  # read in file and generate RPKM and CPM #
  ##########################################
  cat("\nReading in data and checking format ...\n")
  ncol=length(read.table(args[6],sep="\t",nrows=1))
  readcountsRaw=read.table(args[6],header=T,sep="\t",colClasses=c("character",rep("numeric",ncol-1)))
  geneName=readcountsRaw[,1]
  geneLength=readcountsRaw[,2]
  readcountsRaw=readcountsRaw[,-(1:2)]
  rownames(readcountsRaw)=geneName
  
  ncol=length(read.table(args[7],sep="\t",nrows=1))
  colData=read.table(args[7],header=T,sep="\t",colClasses=c("character",rep("factor",ncol-1)))
  # reformat the rownames in design table in case of beginning with a number
  mywarn=FALSE
  for(i in 1:dim(colData)[1]) {
    if(grepl("Processed_",colData[i,1])) {
      colData[i,1]=gsub("Processed_","",colData[i,1])
      mywarn=TRUE
    }
    if(grepl("\\d",substring(colData[i,1],1,1))) {
      colData[i,1]=paste0("X",colData[i,1])
    }
  }
  if(mywarn) { cat("\nDon't include `Processed_` in sample name in the design table next time! Keep the sample name identical!\n")}
  
  rownames(colData)=colData[,1]
  colData=colData[,-1]
  colnames(colData)[1]="condition"
  # check if there is any covariate containing only one level
  mis=0
  for(i in dim(colData)[2]:2){
    if(length(levels(colData[,i]))==1) {
      colData=colData[,-i,drop=F]
      mis=mis+1
    }
  }
  if(mis>1) {
    cat(paste0("\nRemoving ",mis," covariates with only one level ...\n"))
  } else if(mis==1) {
    cat(paste0("\nRemoving ",mis," covariate with only one level ...\n"))
  }
  # check if count and colData are in the same order with respect to samples
  if(sum(rownames(colData)==colnames(readcountsRaw))!=dim(colData)[1]) {
    colData=colData[order(rownames(colData)),]
    readcountsRaw=readcountsRaw[,order(colnames(readcountsRaw))]
  }
  # check if specified contrast interested
  if("contrast"%in%tolower(colnames(colData))) {
    index=which(tolower(colnames(colData))%in%"contrast")
    k0=which(levels(colData$condition)%in%colData[which(colData[,index]=="0")[1],1])
    k1=which(levels(colData$condition)%in%colData[which(colData[,index]=="1")[1],1])
    colData=colData[,-index,drop=F]
  } else {
    k0=1
    head=as.character(lapply(levels(colData$condition), function(x) sub("^([[:alpha:]]*).*","\\1",x)))
    k1=which(head!=head[1])[1]
    if(is.na(k1)) { k1=2}
  }
  
  # generate RPKM file
  rpkm=apply(readcountsRaw,2,function(x) x/sum(x))*10^9/geneLength
  write.csv(rpkm,"RPKM.csv")
  
  # generate CPM file
  cpm=apply(readcountsRaw,2,function(x) x/sum(x)*10^6)
  write.csv(cpm,"CPM.csv")
  cpm=cpm[rowSums(cpm>as.numeric(args[8]))>=min(table(colData[,1])),]
  write.csv(cpm,"filtered_CPM.csv")
  
  ####################################
  # differential expression analysis #
  ####################################
  suppressMessages(dds<-DESeqDataSetFromMatrix(countData=readcountsRaw,colData=colData,design=formula(paste0("~",paste(colnames(colData),sep="",collapse="+")))))
  # filter out lowly expressed genes
  dds=dds[rowSums(fpm(dds,robust=F)>as.numeric(args[8]))>=min(table(colData[,1])),]
  
  cat(paste0("\nDEG analysis using DESeq2\nMaking contrast: ",levels(colData$condition)[k0]," vs ",levels(colData$condition)[k1]),"...\n")
  suppressMessages(dds<-DESeq(dds))
  
  cat("\nWriting results ...\n")
  # Turning off cook's distance testing and independent filtering since we have filtered out lowly expressed genes with cpm
  res=results(dds,contrast=c("condition",levels(colData$condition)[k0],levels(colData$condition)[k1]),alpha=0.05,independentFiltering=F,cooksCutoff=Inf)
  write.csv(res,"DEG_full_list.csv")
  res2=res[res$padj<as.numeric(args[9]) & abs(res$log2FoldChange)>as.numeric(args[10]),]
  write.csv(res2,"DEG_significant_only.csv")
  
  ##################
  # generate plots #
  ##################
  cat("\nTransformation before ploting ...\n")
  # blind=F because we have already performed the function DESeq and then it is not necessary to re-estimate the dispersion values, which will save some time
  rld=rlog(dds,blind=F)
  if(exists("last.warning") && grepl("We recommend instead using the varianceStabilizingTransformation",names(last.warning))) {
    cat("\nNow using the varianceStabilizingTransformation instead ...\n")
    suppressMessages(rld<-vst(dds,blind=F))
  }
  pcaData=plotPCA(rld,intgroup=c("condition"),returnData=TRUE)
  percentVar=round(100*attr(pcaData,"percentVar"))
  
  cat("\nGenerating PCA ...\n")
  suppressMessages(pdf("PCA.pdf",width=9,height=8,paper="special"))
  print(ggplot(pcaData,aes(x=PC1,y=PC2,color=condition))+
          geom_point(size=3)+
          # scale_colour_manual(values=c("#E69F00","#56B4E9"))+
          ggtitle("Principal component analysis")+
          xlab(paste0("PC1: ",percentVar[1],"% variance"))+
          ylab(paste0("PC2: ",percentVar[2],"% variance"))+
          theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.5),
                text=element_text(size=12,family="Helvetica"),
                axis.title=element_text(face="bold"),
                axis.text.x=element_text(size=10,face="bold"),
                axis.text.y=element_text(size=10,face="bold"),
                legend.text=element_text(size=10,face="bold"),
                legend.title=element_text(size=10,face="bold")))
  invisible(dev.off())
  
  suppressMessages(pdf("PCA_with_annotation.pdf",width=9,height=8,paper="special"))
  print(ggplot(pcaData,aes(x=PC1,y=PC2,color=condition))+
          geom_point()+
          ggtitle("Principal component analysis with annotation")+
          geom_text(aes(label=name,fontface=2),hjust=0.5,vjust=-0.8,size=2)+
          xlab(paste0("PC1: ",percentVar[1],"% variance"))+
          ylab(paste0("PC2: ",percentVar[2],"% variance"))+
          theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.5),
                text=element_text(size=12,family="Helvetica"),
                axis.title=element_text(face="bold"),
                axis.text.x=element_text(size=10,face="bold"),
                axis.text.y=element_text(size=10,face="bold"),
                legend.text=element_text(size=10,face="bold"),
                legend.title=element_text(size=10,face="bold")))
  invisible(dev.off())
  
  cat("\nGenerating heatmap ...\n")
  suppressMessages(pdf("heatmap.pdf",width=10,height=8,paper="special",onefile=FALSE))
  pheatmap(cor(assay(rld)),annotation_col=colData,show_colnames=F,main="Heatmap of correlation within samples")
  invisible(dev.off())
  
  cat("\nGenerating volcano plot ...\n")
  category=ifelse(res$padj>as.numeric(args[9])|abs(res$log2FoldChange)<as.numeric(args[10]),"Regularly expressed",ifelse(res$log2FoldChange>0,paste0("Highly expressed in ",levels(colData[,1])[k0]),paste0("Highly expressed in ",levels(colData[,1])[k1])))
  index=as.factor(category)
  if(length(levels(index))==3) {
    index=factor(category,levels=c("Regularly expressed",paste0("Highly expressed in ",levels(colData[,1])[k0]),paste0("Highly expressed in ",levels(colData[,1])[k1])),ordered=T)
  }
  plot=data.frame(logFC=res$log2FoldChange,log10p=-log10(res$padj),index=index)
  suppressMessages(pdf("volcano_plot.pdf",width=10,height=8,paper="special"))
  print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
          geom_point(aes(size=index))+
          ggtitle(paste0("Volcano plot of ",levels(colData$condition)[k0]," vs ",levels(colData$condition)[k1]))+
          theme_bw()+theme_classic()+
          scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
          scale_size_manual(values=c(1,2,2))+
          scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
          geom_vline(xintercept=as.numeric(args[10]),colour="black",linetype=2)+
          geom_vline(xintercept=-as.numeric(args[10]),colour="black",linetype=2)+
          geom_hline(yintercept=-log10(as.numeric(args[9])),colour="black",linetype=2)+
          theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.5),
                text=element_text(size=12,family="Helvetica"),
                axis.title=element_text(face="bold"),
                axis.text.x=element_text(size=10,face="bold"),
                axis.text.y=element_text(size=10,face="bold"),
                legend.text=element_text(size=10,face="bold"),
                legend.title=element_text(size=10,face="bold")))
  invisible(dev.off())
  
  cat("\nGenerating MA plot ...\n")
  suppressMessages(pdf("MA_plot.pdf",width=8,height=8,paper="special"))
  plotMA(dds,ylim=c(-3,3),main="DESeq2 MA plot")
  invisible(dev.off())
  
  #################
  # Gene Ontology #
  #################
  if(args[11]=="TRUE" & dim(res2)[1]>0) {
    cat("\nPerforming Gene Ontology analysis ...\n")
    xegdb=paste0("org.",args[12],".eg.db")
    xegGO=paste0("org.",args[12],".egGO")
    # extract Entrez IDs
    suppressMessages(E_IDs<-select(get(xegdb),rownames(res2),c("ENTREZID","GENENAME"),"ALIAS"))
    Entrez=unique(E_IDs[!is.na(E_IDs$ENTREZID),2])
    universe=mappedkeys(get(xegGO))
    
    if(sum(Entrez%in%universe)>0) {
      Entrez=Entrez[Entrez%in%universe]
      paramsBP=new('GOHyperGParams',geneIds=Entrez,universeGeneIds=universe,ontology='BP',pvalueCutoff=as.numeric(args[13]),conditional=F,testDirection='over',annotation=xegdb)
      hgOverBP=hyperGTest(paramsBP)
      resultBP=summary(hgOverBP)
      write.csv(resultBP,file="GO_BP.csv")
      htmlReport(hgOverBP,file="GO_BP.html",digits=5)
      
      paramsCC=new('GOHyperGParams',geneIds=Entrez,universeGeneIds=universe,ontology='CC',pvalueCutoff=as.numeric(args[13]),conditional=F,testDirection='over',annotation=xegdb)
      hgOverCC=hyperGTest(paramsCC)
      resultCC=summary(hgOverCC)
      write.csv(resultCC,file="GO_CC.csv")
      htmlReport(hgOverCC,file="GO_CC.html",digits=5)
      
      paramsMF=new('GOHyperGParams',geneIds=Entrez,universeGeneIds=universe,ontology='MF',pvalueCutoff=as.numeric(args[13]),conditional=F,testDirection='over',annotation=xegdb)
      hgOverMF=hyperGTest(paramsMF)
      resultMF=summary(hgOverMF)
      write.csv(resultMF,file="GO_MF.csv")
      htmlReport(hgOverMF,file="GO_MF.html",digits=5)
      
      paramsKEGG=new('KEGGHyperGParams',geneIds=Entrez,universeGeneIds=universe,pvalueCutoff=as.numeric(args[13]),testDirection='over',annotation=xegdb)
      hgOverKEGG=hyperGTest(paramsKEGG)
      resultKEGG=summary(hgOverKEGG)
      write.csv(resultKEGG,file="KEGG.csv")
    } else {
      cat("\nGene IDs of DEGs cannot be identified in the GOgeneIds!!!\nIs the SPECIES correctly specified?")
    }
  } else if(dim(res2)[1]==0) {
    cat("\nNo significant DEGs! Skip GO analysis ...\n")
  } else {
    cat("\nGO analysis is not specified. Skip GO analysis ...\n")
  }
  cat("\nDEG analysis finished!\n\n")
}
