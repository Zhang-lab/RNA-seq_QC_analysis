#!/usr/bin/env Rscript
args=commandArgs()

name=args[6]
genome=args[7]
beta=args[8]
rtime=args[9]
image_id=args[10]

# TXT report##############################################################################################
report=list(Library=name)

dedup=read.table(paste("dedup_percentage_",name,".result",sep=""),header=T,sep="\t")
if(dim(dedup)[1]>1) {
  dedup=data.frame(filename=strsplit(as.character(dedup[1,1]),'_1',fixed=T)[[1]][1],ratio=100-(dedup[1,2]+dedup[2,2])/2,class='Sample')
  data_type="Paired-end data"
} else {
  dedup=data.frame(filename=strsplit(as.character(dedup[1,1]),'_1',fixed=T)[[1]][1],ratio=100-dedup[1,2],class='Sample')
  data_type="Single-end data"
}

report=append(report,list(Pipeline.version=args[8],Genome=genome,Data.type=data_type))
##############################################################################################
qc=read.table('RNA-seq_qc_collection.txt',header=T)
sample=as.numeric(qc[,c(2,3,4,5,6)])

library=data.frame(Sample=sample)
rownames(library)=c("Total reads","Written reads by cutadapt","Reads aligned 0 time","Reads aligned exactly 1 time","Reads aligned greater than 1 time")
report=append(report,list(Library.size=library))
##############################################################################################
sample=c(paste(round(dedup[1,2],2),"%",sep=""),paste(round(qc[1,6]*100.0/sum(qc[1,4:6]),2),"%",sep=""))
library=data.frame(Sample=sample)
rownames(library)=c("Before alignment library duplicates percentage","After alignment PCR duplicates percentage")
report=append(report,list(Library.complexity=library))
##############################################################################################
sample=c(paste(round(qc[,7]*100,2),"%",sep=""),as.numeric(qc[,8]),as.numeric(qc[,9]))
library=data.frame(Sample=sample)
rownames(library)=c("Rate of reads with gene feature","Area under curve of gene body coverage","Expressed genes with CPM>1")
report=append(report,list(Feature.count=library))
##############################################################################################
read=read.table(paste("step3.6",name,"reads_distri_table.txt",sep="_"),header=T)
sample=read[,3]
library=data.frame(sample)
colnames(library)="Tag counts"
rownames(library)=as.character(read[,1])
report=append(report,list(Read.distribution=library))
##############################################################################################
gc=read.table(paste("step3.5_",name,".GC.txt",sep=""),header=T)
sample=c(sum(gc[which(gc[,1]<25),2]),sum(gc[which(gc[,1]>=25&gc[,1]<50),2]),sum(gc[which(gc[,1]>=50&gc[,1]<75),2]),sum(gc[which(gc[,1]>=75&gc[,1]<=100),2]))
library=data.frame(sample)
colnames(library)="Read counts"
rownames(library)=c("GC in 0-25","GC in 25-50","GC in 50-75","GC in 75-100")
report=append(report,list(GC.content=library))
##############################################################################################
capture.output(print(report),file=paste(name,"report.txt",sep='_'))

# JSON report##############################################################################################
if(is.na(image_id)) {
  part1=data.frame(name,genome,data_type,beta,"MD5ToBeChange",rtime)
  colnames(part1)=c("file_name","genome","read_type","pipe_version","bash_script_MD5","running_time")
  file=list(`data_information`=part1)
} else {
  part1=data.frame(name,genome,data_type,beta,image_id,"MD5ToBeChange",rtime)
  colnames(part1)=c("file_name","genome","read_type","pipe_version","Docker_image_id","bash_script_MD5","running_time")
  file=list(`data_information`=part1)
}

part2=data.frame("cutadapt","1.16",as.numeric(qc[,3]),"FastQC","0.11.5")
colnames(part2)=c("program1","program1_version","written_reads_by_cutadapt","program2","program2_version")
file=append(file,list(`pre_alignment_stats`=part2))

part3=data.frame("STAR","2.5.4b","--outWigType bedGraph --outWigNorm RPM --outSAMtype BAM SortedByCoordinate",as.numeric(qc[,4]),as.numeric(qc[,5]),as.numeric(qc[,6]))
colnames(part3)=c("alignment_program","alignment_program_version","alignment_program_parameters","reads_aligned_0_time","reads_aligned_exactly_1_time","reads_aligned_greater_than_1_time")
file=append(file,list(`mapping_stats`=part3))

part9=data.frame(round(dedup[1,2],2)/100,round(qc[1,6]*100.0/sum(qc[1,4:6]),2)/100)
colnames(part9)=c("before_alignment_library_duplicates_percentage","after_alignment_PCR_duplicates_percentage")
file=append(file,list(`library_complexity`=part9))

gene_type=read.table(paste0("step4.2_gene_type_count_",name,".txt"),header=T)
gene_type[,1]=paste0("@",gene_type[,1],"@")
gene_type=paste(gene_type[,1],gene_type[,3],sep=": ")
gene_type=paste0("!",paste(gene_type,sep="",collapse=", "),"!")

part4=data.frame("featureCounts","1.5.1","-p -O -g <gene_name>/<gene_type>",as.numeric(round(qc[,7],4)),as.numeric(qc[,8]),as.numeric(qc[,9]),gene_type)
colnames(part4)=c("software","version","count_parameter","rate_of_reads_with_gene_feature","area_under_curve_of_gene_body_coverage","expressed_genes_with_cpm_larger_than_1","gene_type_fragment_count")
file=append(file,list(`feature_counting`=part4))

yield=read.table(paste0("yield_",name,".result"),sep='\t',header=T)
yield=yield[yield$TOTAL_READS<=1e8,]

part10=data.frame(paste("?",paste(yield[,1],sep="",collapse=","),"?",sep=""),paste("?",paste(yield[,2],sep="",collapse=","),"?",sep=""))
colnames(part10)=c("total_reads","expected_distinction")
file=append(file,list(`yield_distribution`=part10))

gene_cov=read.table(paste0("step3.1_",name,".geneBodyCoverage.txt"))
part11=data.frame(percentile=paste0("?",paste(gene_cov[1,-1],sep="",collapse=","),"?"),coverage=paste0("?",paste(gene_cov[2,-1],sep="",collapse=","),"?"))
file=append(file,list(`gene_body_covergae`=part11))

gcon=data.frame(c("GC_in_0-25","GC_in_25-50","GC_in_50-75","GC_in_75-100"),c(sum(gc[which(gc[,1]<25),2]),sum(gc[which(gc[,1]>=25&gc[,1]<50),2]),sum(gc[which(gc[,1]>=50&gc[,1]<75),2]),sum(gc[which(gc[,1]>=75&gc[,1]<=100),2])))
gcon[,1]=paste0("@",gcon[,1],"@")
gcon=paste(gcon[,1],gcon[,2],sep=": ")
gcon=paste0("!",paste(gcon,sep="",collapse=", "),"!")

rd=data.frame(as.character(read[,1]),as.numeric(read[,3]))
rd[,1]=paste0("@",rd[,1],"@")
rd=paste(rd[,1],rd[,2],sep=": ")
rd=paste0("!",paste(rd,sep="",collapse=", "),"!")

junction=read.table(paste0("step3.8_",name,".result"),header=F)
junction[,1]=paste0("@",junction[,1],"@")
event=paste0("!",paste(paste(junction[,1],junction[,2],sep=": "),sep="",collapse=", "),"!")
junction=paste0("!",paste(paste(junction[,1],junction[,3],sep=": "),sep="",collapse=", "),"!")

part6=data.frame(GC_content=gcon,reads_distribution=rd,splice_junction=junction,splice_events=event)
file=append(file,list(`RseQC_report`=part6))

saturate=read.table(paste0("step4.3_CPM_saturation_",name,".txt"),header=T)
part5=data.frame(t(apply(saturate,2,function(x) paste0("?",paste(x,sep="",collapse=","),"?"))))
colnames(part5)=c("sequence_depth","cpm_from_1_to_10","cpm_from_10_to_50","cpm_greater_than_50","total genes")
file=append(file,list(`saturation`=part5))

library(jsonlite)
capture.output(toJSON(file,pretty=T),file=paste(name,"report.json",sep='_'))

