#!/bin/bash  
species=$1

# get pipe path, though readlink/realpath can do it, some version doesn't have that
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
 
fastq_dump_tool='fastq-dump.2.8.2'
adapter_1="AGATCGGAAGAGCACACGTCTGAAC"
adapter_2="AGATCGGAAGAGCGTCGTGTAGGGA"
preseq="preseq"
RSeQC_script="/home/shaopengliu/test_rna/scripts"
ref_gene_bed="/home/shaopengliu/test_rna/mm10_GENCODE_VM11_basic.bed"
gencode_region='/home/shaopengliu/resources/mm10/mm10_gencode_vM15_gtf_region.bed'


# genome specific resources:
if [[ $species == mm10 ]]; 
	then
	#star_ref='/home/Resource/Genome/mm10/STAR_index_mm10.gencode.vM9'  
	annotation_file="/home/shaopengliu/resources/mm10/gencode.vM15.annotation.gtf" 
	chrom_size="/home/Resource/Genome/mm10/mm10.chrom.sizes"
#elif [[ $species == mm9 ]];
#	then
#
#elif [[ $species == hg38 ]];
#	then
#
#elif [[ $species == hg19 ]];
#	then
#
#elif [[ $species == danRer10 ]];
#	then
#
fi

# other tools note (their names are stable, so I use it directly)
# 1, cutadapt v1.12
# 2, fastqc v0.11.5
# 3, multiqc v1.2
# 4, preseq 2.0.0
# 5, STAR 2.5.3a
