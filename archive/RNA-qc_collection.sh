#!/bin/bash
name=$1
types=$2
chrom_size='/home/Resource/Genome/mm10/mm10.chrom.sizes'
gencode_region='/home/shaopengliu/resources/mm10/mm10_gencode_vM15_gtf_region.bed'
annotation_file='/home/shaopengliu/resources/mm10/gencode.vM15.annotation.gtf'
threads=24
RSeQC_script='/home/shaopengliu/test_rna/scripts'
ref_gene_bed='/home/shaopengliu/test_rna/mm10_GENCODE_VM11_basic.bed'
adapter_1="AGATCGGAAGAGCACACGTCTGAAC"
adapter_2="AGATCGGAAGAGCGTCGTGTAGGGA"

if [[ $types == PE ]];
	then
	echo 'trimming RNA-seq PE reads by cutadapt'
	cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=36  -o 'Trimed_'$name'_1.fastq' -p  'Trimed_'$name'_2.fastq'  *fastq.gz  > 'step1.1_'$name'_cutadapt_PE.trimlog'  
	temp=`grep "Total read pairs processed:" step1.1_*trimlog | awk '{print $5}'`  
	raw_reads=`echo ${temp//,}`  
elif [[ $types == SE ]];
	then
	echo 'trimming RNA-seq SE reads by cutadapt'
	cutadapt -a $adatper_1 --quality-cutoff=15,10 --minimum-length=36  -o  'Trimed_'$name'.fastq' *fastq.gz > 'step1.1_'$name'_cutadapt_SE.trimlog'  
	temp=`grep "Total reads processed:" step1.1_*trimlog | awk '{print $4}'`
	raw_reads=`echo ${temp//,}` 
fi


# 1,  total reads from cutadapt trimlog:
if [[ $types == PE ]];
	then
	temp=`grep "Total read pairs processed:" step1.1_*trimlog | awk '{print $5}'`  
	raw_reads=`echo ${temp//,}`  
	temp2=`grep "Pairs that were too short" step1.1_*trimlog | awk '{print $6}'`
	removed_reads=`echo ${temp2//,}`
elif [[ $types == SE ]];
	then
	temp=`grep "Total reads processed:" step1.1_*trimlog | awk '{print $4}'`
	raw_reads=`echo ${temp//,}` 
	temp2=`grep "Reads that were too short" step1.1_*trimlog | awk '{print $6}'`
	removed_reads=`echo ${temp2//,}`
fi

#2,library duplicate
# note: it seems that the HISAT2 would record the multiple match in bam, so the total number in methylQA output is larger than the original raw reads
# regarding the 
#2.1, from fastqc:
# see "dedup percentage in result file"

#2.2, from mapping: Hisat2 output
if [[ $types == PE ]];
then
hisat2 -q -x /home/shaopengliu/resources/mm10/hisat2_index/ht2_idx_mm10 -1 'Trimed_'$name'_1.fastq'  -2 'Trimed_'$name'_2.fastq'  2> step2.1_hisat2_summary.txt | samtools view -bS - | samtools sort - -o step2.1_hisat2_sorted_$name'.bam'
echo "step2.1, mapping as PE" >> pipe_processing.log
else
hisat2 -q -x /home/shaopengliu/resources/mm10/hisat2_index/ht2_idx_mm10 -U 'Trimed_'$name'.fastq' 2> step2.1_hisat2_summary.txt | samtools view -bS - | samtools sort - -o step2.1_hisat2_sorted_$name'.bam'
echo "step2.1, mapping as SE >> pipe_processing.log"
fi

# step2.1_hisat2_summary.txt
align_0_time_num=`grep "0 times" step2.1_hisat2_summary.txt | head -1 | awk '{print $1}'`
exact_1_time_num=`grep "exactly 1 time" step2.1_hisat2_summary.txt | head -1 | awk '{print $1}'`
gt_1_time_num=`grep ">1 times" step2.1_hisat2_summary.txt | head -1 | awk '{print $1}'`


#4,  genebody coverage area under curve:
rm step3.1*
python2.7 $RSeQC_script'/geneBody_coverage.py' -i step2.1_hisat2_sorted_$name'.bam' -r $ref_gene_bed -o 'step3.1_'$name 
auc=`head -1 step3.1_*.r | awk -F "[(]" '{print $2}' | sed 's/)//' | awk -F "[,]" '{for(i=1;i<=NF;i++) t+=$i; print t}'`

#5, gene type count
rm step4.*
if [[ $types == PE ]];
	then
	echo 'featureCounts on PE data'
	featureCounts -a $annotation_file -p -T $threads -O  \
	 -o step4.1_gene_name_fc_$name  -g gene_name step2.1_hisat2_sorted_$name'.bam'
	featureCounts -a $annotation_file -p -T $threads -O  \
	 -o step4.1_gene_type_fc_$name -g gene_type step2.1_hisat2_sorted_$name'.bam'
elif [[ $types == SE ]];
	then
	echo 'featureCounts on SE data'
	featureCounts -a $annotation_file -T 24 -O  \
	 -o step4.1_gene_name_fc_$name  -g gene_name step2.1_hisat2_sorted_$name'.bam'
	featureCounts -a $annotation_file  -T 24 -O  \
	 -o step4.1_gene_type_fc_$name -g gene_type step2.1_hisat2_sorted_$name'.bam'
fi

# 4.2, FC collection
echo -e "gene_name\tlength\tfragment_count" > 'step4.2_gene_name_count_'$name'.txt'
tail -n +3 step4.1_gene_name_fc_$name | awk '{print $1"\t"$6"\t"$7}' | sort -k3,3rn >> 'step4.2_gene_name_count_'$name'.txt'

echo -e "gene_type\tlength\tfragment_count" > 'step4.2_gene_type_count_'$name'.txt'
tail -n +3 step4.1_gene_type_fc_$name | awk '{print $1"\t"$6"\t"$7}' | sort -k3,3rn >> 'step4.2_gene_type_count_'$name'.txt'

cp step4.2_gene_type_count*txt ./data_collection_$name

#3, reads with gene feature:
total_count=`awk '{s+=$2}END{print s}' step4.1_gene_name_fc*summary`
reads_with_feature=`grep Assigned step4.1_gene_name_fc*summary | awk '{print $2}'`
rates_with_feature=`python -c "print(1.0*$reads_with_feature/$total_count)"`


#6, reads distribution in Genome
cp step3.6_*read_distribution.txt ./data_collection_$name

#7, GC content
cp step3.5_*.GC.xls  ./data_collection_$name

#8, expression matrix and CPM
cp step4.2_gene_name*txt ./data_collection_$name
total_count=`awk '{s+=$3}END{print s}' step4.2_gene_name*txt`
thres=`python -c "print(1.0* $total_count/1000000)"`
exp_gene_num=`awk -v thres=$thres '$3>thres' step4.2_gene_name*txt | wc -l`

#summarize together:
echo -e "name\traw_reads\tremoved_reads\talign_0_time\texactly_1_time\tgt_1_time\treads_ratio_with_gene_feature\tauc_genebody\texp_gene_num" > RNA-seq_qc_collection.txt
echo -e "$name\t$raw_reads\t$removed_reads\t$align_0_time_num\t$exact_1_time_num\t$gt_1_time_num\t$rates_with_feature\t$auc\t$exp_gene_num" >> RNA-seq_qc_collection.txt
mv RNA-seq_qc_collection.txt ./data_collection_$name

