#!/usr/bin/env Rscript
args=commandArgs()

name=args[6]
genome=args[7]
beta=args[8]
rtime=args[9]
image_id=args[10]

# TXT report##############################################################################################
report=list(Library=name)

dedup=read.table(paste("step1.3_dedup_percentage_",name,".result",sep=""),header=T,sep="\t")
if(dim(dedup)[1]>1) {
  dedup=data.frame(filename=strsplit(as.character(dedup[1,1]),'_1',fixed=T)[[1]][1],ratio=100-(dedup[1,2]+dedup[2,2])/2,class='Sample')
  data_type="Paired-end data"
} else {
  dedup=data.frame(filename=strsplit(as.character(dedup[1,1]),'_1',fixed=T)[[1]][1],ratio=100-dedup[1,2],class='Sample')
  data_type="Single-end data"
}

report=append(report,list(Pipeline.version=args[8],Genome=genome,Data.type=data_type))
##############################################################################################
qc=read.table(paste0("QC_summary_",name,".txt"),header=T)
sample=as.numeric(qc[1,c(2,3,3,4,7,8,10)])
sample[3]=paste0(round(sample[3]/sample[1],4)*100,"%")

library=data.frame(Sample=sample)
rownames(library)=c("Total reads","Written reads by cutadapt","Written reads percentage","Reads average length","Reads aligned 0 time","Reads aligned exactly 1 time","Reads aligned greater than 1 time")
report=append(report,list(Library.size=library))
##############################################################################################
sample=c(paste0(round(qc[1,5]*100,2),"%"),paste0(round(qc[1,6]*100,2),"%"))
library=data.frame(Sample=sample)
rownames(library)=c("Before alignment library duplicates percentage","After alignment duplicates percentage")
report=append(report,list(Library.complexity=library))
##############################################################################################
sample=c(paste(round(qc[1,13]*100,2),"%",sep=""),as.numeric(qc[1,11]),as.numeric(qc[1,14]),as.numeric(qc[1,15]))
library=data.frame(Sample=sample)
rownames(library)=c("Ratio of reads with gene feature","Area under curve of gene body coverage","All genes number","Expressed genes with CPM>1")
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
sample=c(as.character(median(gc[,1])),sum(gc[which(gc[,1]<25),2]),sum(gc[which(gc[,1]>=25&gc[,1]<50),2]),sum(gc[which(gc[,1]>=50&gc[,1]<75),2]),sum(gc[which(gc[,1]>=75&gc[,1]<=100),2]))
library=data.frame(Sample=sample)
rownames(library)=c("Median GC content","Read counts of GC in 0-25","Read counts of GC in 25-50","Read counts of GC in 50-75","Read counts of GC in 75-100")
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

part2=data.frame("cutadapt","1.16",as.numeric(qc[1,2]),as.numeric(qc[1,3]),round(qc[1,3]/qc[1,2],4),"FastQC","0.11.5")
colnames(part2)=c("program1","program1_version","total_reads","written_reads_by_cutadapt","written_reads_ratio","program2","program2_version")
file=append(file,list(`pre_alignment_stats`=part2))

splice=read.table(paste0("step2.1_Star_Splice_",name,".txt"),sep="\t",colClasses=c("character","numeric"))
splice[,1]=unlist(lapply(splice[,1],function(x) unlist(strsplit(x,":"))[2]))
splice[,1]=paste0("@",splice[,1],"@")
splice=paste(splice[,1],splice[,2],sep=": ")
splice=paste0("!",paste(splice,collapse=", "),"!")

part3=data.frame("STAR","2.5.4b","--outWigType bedGraph --outWigNorm RPM --outSAMtype BAM SortedByCoordinate --outWigStrand Stranded --quantMode TranscriptomeSAM",as.numeric(qc[1,4]),as.numeric(qc[1,7]),as.numeric(qc[1,8]),as.numeric(qc[1,9]),as.numeric(qc[1,10]),splice)
colnames(part3)=c("alignment_program","alignment_program_version","alignment_program_parameters","average_reads_length","reads_aligned_0_time","reads_aligned_exactly_1_time","uniquely_mapped_ratio","reads_aligned_greater_than_1_time","number_of_splice")
file=append(file,list(`mapping_stats`=part3))

part9=data.frame(round(qc[1,5:6],4))
colnames(part9)=c("before_alignment_library_duplicates_ratio by FastQC v0.11.5","after_alignment_duplicates_ratio by Picard v1.3.2")
file=append(file,list(`library_complexity`=part9))

gene_type=read.table(paste0("step4.2_gene_type_count_",name,".txt"),header=T)
gene_type[,1]=paste0("@",gene_type[,1],"@")
gene_type=paste(gene_type[,1],gene_type[,3],sep=": ")
gene_type=paste0("!",paste(gene_type,sep="",collapse=", "),"!")

saturate=read.table(paste0("step4.3_CPM_saturation_",name,".txt"),header=T)
ex.dis=paste0("!",paste(c("@between_1_and_10@","@between_10_and_50@","@larger_than_50@"),as.numeric(saturate[nrow(saturate),2:4]),sep=": ",collapse=", "),"!")

type_distri=read.table(paste0("step4.2_temp_store_",name,".txt"), header=F)
rownames(type_distri)<-type_distri[, 1]

part4=data.frame("featureCounts","1.5.1","gtf_ToBeChange","-p -g <gene_name>/<gene_type> -Q 10 -s <0,1,2>",as.numeric(round(qc[1,12],4)),as.numeric(round(qc[1,13],4)),as.numeric(qc[1,11]),  type_distri["Mt_rRNA", 2],type_distri["Mt_tRNA", 2],type_distri["ribozyme", 2], type_distri["rRNA", 2]    ,as.numeric(qc[1,14]),as.numeric(qc[1,15]),ex.dis,gene_type)
colnames(part4)=c("software","version","annotation_file_version","count_parameter","number_of_uniquely_mapped_fragments_with_gene_feature","ratio_of_uniquely_mapped_fragments_with_gene_feature","area_under_curve_of_gene_body_coverage",  "reads_in_Mt_rRNA", "reads_in_Mt_tRNA", "reads_in_ribozyme", "reads_in_rRNA"  ,"detected_genes_number","detected_genes_with_cpm_larger_than_1","detected_genes_cpm_distribution","gene_type_fragment_count")
file=append(file,list(`feature_counting`=part4))

if (file.exists(paste0("step4.2_temp_globin_",name,".txt"))) {
  globin_count=read.table(paste0("step4.2_temp_globin_",name,".txt"), header=F)
  part42<-data.frame(t(data.frame(globin_count$V2, row.names=globin_count$V1)), row.names=NULL)
  file=append(file, list(`globin_genes`=part42))
}

yield=read.table(paste0("step2.3_yield_",name,".result"),sep='\t',header=T)
yield=yield[yield$TOTAL_READS<=1e8,]

part10=data.frame(paste("?",paste(yield[,1],sep="",collapse=","),"?",sep=""),paste("?",paste(yield[,2],sep="",collapse=","),"?",sep=""))
colnames(part10)=c("total_reads","expected_distinction")
file=append(file,list(`yield_distribution`=part10))

gene_cov=read.table(paste0("step3.1_",name,".geneBodyCoverage.txt"))
part11=data.frame(percentile=paste0("?",paste(gene_cov[1,-1],sep="",collapse=","),"?"),coverage=paste0("?",paste(gene_cov[2,-1],sep="",collapse=","),"?"))
file=append(file,list(`gene_body_covergae`=part11))

gc.index=c()
gc.count=c()
for (i in seq(0,96,2)) {
  gc.index=c(gc.index,paste0("GC_in_",i,"-",i+2))
  gc.count=c(gc.count,sum(gc[which(gc[,1]>=i & gc[,1]<(i+2)),2]))
}
gc.index=c(gc.index,"GC_in_98-100")
gc.count=c(gc.count,sum(gc[which(gc[,1]>=98 & gc[,1]<=100),2]))
gcon=data.frame(gc.index,gc.count)
gcon[,1]=paste0("@",gcon[,1],"@")
gcon=paste(gcon[,1],gcon[,2],sep=": ")
gcon=paste0("!",paste(gcon,sep="",collapse=", "),"!")

rd=data.frame(as.character(read[,1]),as.numeric(read[,3]))
rd[,1]=paste0("@",rd[,1],"@")
rd=paste(rd[,1],rd[,2],sep=": ")
rd=paste0("!",paste(rd,sep="",collapse=", "),"!")

junction=read.table(paste0("step3.4_",name,".result"),header=F)
junction[,1]=paste0("@",junction[,1],"@")
event=paste0("!",paste(paste(junction[,1],junction[,3],sep=": "),sep="",collapse=", "),"!")
junction=paste0("!",paste(paste(junction[,1],junction[,2],sep=": "),sep="",collapse=", "),"!")

part6=data.frame(median_GC_content=median(gc[,1]),GC_content=gcon,reads_distribution=rd,splice_junction=junction,splice_events=event)
file=append(file,list(`RseQC_report`=part6))

part5=data.frame(t(apply(saturate,2,function(x) paste0("?",paste(x,sep="",collapse=","),"?"))))
colnames(part5)=c("gene_count_subsampling","cpm_from_1_to_10","cpm_from_10_to_50","cpm_greater_than_50","total_genes")
file=append(file,list(`saturation`=part5))

library(jsonlite)
capture.output(toJSON(file,pretty=T),file=paste(name,"report.json",sep='_'))


