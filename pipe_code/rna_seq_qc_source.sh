#!/bin/bash  
species=$1

# get pipe path, though readlink/realpath can do it, some version doesn't have that
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
 
# executables
cutadapt="/usr/bin/cutadapt"
RSeQC_script="/home/shaopengliu/pipe_script/github/RNA-seq_QC_analysis/pipe_code/scripts"
RSeQC_python="python2.7"  #for simg consistency
fastq_dump_tool='/usr/bin/fastq-dump.2.8.1-3'
preseq="/usr/bin/preseq"
fastqc="/usr/bin/fastqc"
STAR="/usr/bin/STAR"
samtools="/usr/local/bin/samtools"
featureCounts="/usr/bin/featureCounts"
rsem="/home/shaopengliu/tools/anaconda3/bin/rsem-calculate-expression"
multiqc="/usr/local/bin/multiqc"

# genome specific resources:
if [[ $species == mm10 ]]; 
    then
    star_ref='/home/Resource/Genome/mm10/STAR_index_mm10_v2.5.4b.gencode.vM15'
    echo "the specified genome is mm10"  
    annotation_file="/home/shaopengliu/resources/mm10/gencode.vM15.annotation.gtf" 
    chrom_size="/home/Resource/Genome/mm10/mm10.chrom.sizes"
    ref_gene_bed="/home/shaopengliu/resources/mm10/mm10_GENCODE_VM11_basic.bed"
    sub_ref_gene="/home/shaopengliu/resources/mm10/mm10_GENCODE_VM11_basic_subsample.bed"
    rsem_ref="/home/shaopengliu/resources/mm10/rsem_ref/mm10_rsem_ref"
elif [[ $species == hg38 ]];
    then
    echo "the specified genome is hg38"
    star_ref='/home/Resource/Genome/hg38/STAR_index_gencode.v27.annotation'
    annotation_file="/home/shaopengliu/resources/hg38/gencode.v27.annotation.gtf" 
    chrom_size="/home/Resource/Genome/hg38/hg38.25_chromsome.sizes"
    ref_gene_bed="/home/shaopengliu/resources/hg38/hg38_RefSeq.bed"
    sub_ref_gene="/home/shaopengliu/resources/hg38/hg38_RefSeq_subsample_topexp.bed"
    rsem_ref="/home/shaopengliu/resources/hg38/rsem_ref/hg38_rsem_ref"
elif [[ $species == hg19 ]];
    then
    echo "the specified genome is hg19"
    star_ref="/home/Resource/Genome/hg19/STAR_index_hg19_gencodeV24"
    annotation_file="/home/shaopengliu/resources/hg19/gencode.v24lift37.annotation.gtf"
    chrom_size="/home/Resource/Genome/hg19/hg19_chromosome.size"
    ref_gene_bed="/home/shaopengliu/resources/hg19/hg19_RefSeq_RSeQC.bed"
    rsem_ref="/home/shaopengliu/resources/hg19/rsem_ref/hg19_rsem_ref"
elif [[ $species == danRer10 ]];
    then
    echo "the specified genome is danRer10"
    star_ref="/home/Resource/Genome/danRer10/STAR_index_GRCz10_83"
    annotation_file="/home/shaopengliu/resources/danRer10/Danio_rerio.GRCz10.83.gtf"
    chrom_size="/home/shaopengliu/resources/danRer10/danRer10.chrom.sizes"
    ref_gene_bed="/home/shaopengliu/resources/danRer10/Refgene_danRer10_RSeQC.bed"
    rsem_ref="/home/shaopengliu/resources/danRer10/rsem_ref/danRer10_rsem_ref"
fi

if [ -z "$sub_ref_gene" ]
    then
    sub_ref_gene=$ref_gene_bed
fi



