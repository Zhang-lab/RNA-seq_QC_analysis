#!/bin/bash
# pipe usage:
# user@domain: path_to_pipe/pipe.sh -g <hg38/mm10/danRer10>  -r <PE/SE>  -o file1  -p file2

# pipe start
###################################################################################################
# Preparation:
date

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 

# read parameters
while getopts m:t:g:o:p:r:  opts
do case "$opts" in
m) marker="$OPTARG";;	# default 'unmarked'
t) threads="$OPTARG";;	# default 24
g) species="$OPTARG";;	# hg19, hg38, mm9, mm10, danRer10
o) R1="$OPTARG";;    # PE read 1, or the SE file, or the sra file
p) R2="$OPTARG";;    # PE read 2. 
r) types="$OPTARG";;	# PE or SE;
h) echo "usage:  path-to-pipe/pipe.sh  -g <hg38/hg19/mm10/mm9/danRer10/personal>  -r <PE/SE> -o read_file1  -p read_file2 (if necessary)"
exit;;
[?]) "Usage: ./pipe.sh  -g <hg38/hg19/mm10/mm9/danRer10> -o <read_file1>  -r <PE/SE>";;
esac
done

if [ -z "$threads" ]
then
threads=24
fi

if [ -z "$marker" ]
then
echo "you didn't specify the marker, are you sure to keep the default unmarked"
marker='unmarked'
fi

source $pipe_path'/rna_seq_qc_source.sh' $species

if [[ $R1 == *.sra ]]
	then name=`echo ${R1%.sra}`
	echo "this is sra file, $fastq_dump_tool would be used......"
	$fastq_dump_tool $R1 --split-3
	raw1=$name'_1.fastq'
	raw2=$name'_2.fastq'
elif [[ $R1 == *.fastq* ]] && [[ $types == PE  ]]
	then
	name=`echo ${R1%.fastq*}`
	raw1=$R1
	raw2=$R2
elif [[ $R1 == *.fastq* ]] && [[ $types == SE ]]
	then
	name=`echo ${R1%.fastq*}`
	raw1=$R1
else
	echo "please use fastq(or fastq.gz) file or sra file......"
	exit
fi


mkdir 'Processed_'$name
mv $R1    ./'Processed_'$name/
mv pesudo_bl.txt  ./'Processed_'$name/  2> /dev/null
mv $raw1  ./'Processed_'$name/  2> /dev/null
mv $raw2  ./'Processed_'$name/  2> /dev/null
cd ./'Processed_'$name/
mkdir 'data_collection_'$name
touch pipe_processing.log

# refine chrom_size file (remove random and Unknown record)
awk  '{ if ((length($1) < 6) && (length($1) > 1))  print $0}' OFS='\t' $chrom_size  > refined_chrom_size.txt
chrom_size=`pwd`"/refined_chrom_size.txt"

# start record
date >> pipe_processing.log
echo "Target file is $R1 $R2" >> pipe_processing.log
echo "Specified species is $species" >> pipe_processing.log
echo "types of reads is $types" >> pipe_processing.log
echo " " >> pipe_processing.log

###################################################################################################
# Step 1, Trim ATAC-seq adapters and QC on seq file
# 1.1 Trim by cutadapt
if [[ $types == PE ]];
	then
	echo 'trimming RNA-seq PE reads by cutadapt'
	cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=36  -o 'Trimed_'$name'_1.fastq' -p  'Trimed_'$name'_2.fastq'  $raw1 $raw2  > 'step1.1_'$name'_cutadapt_PE.trimlog'  
	temp=`grep "Total read pairs processed:" step1.1_*trimlog | awk '{print $5}'`  
	raw_reads=`echo ${temp//,}`  
	temp2=`grep "Pairs that were too short" step1.1_*trimlog | awk '{print $6}'`
	removed_reads=`echo ${temp2//,}`
elif [[ $types == SE ]];
	then
	echo 'trimming RNA-seq SE reads by cutadapt'
	cutadapt -a $adatper_1 --quality-cutoff=15,10 --minimum-length=36  -o  'Trimed_'$name'.fastq' $raw1 > 'step1.1_'$name'_cutadapt_SE.trimlog'  
	temp=`grep "Total reads processed:" step1.1_*trimlog | awk '{print $4}'`
	raw_reads=`echo ${temp//,}` 
	temp2=`grep "Reads that were too short" step1.1_*trimlog | awk '{print $6}'`
	removed_reads=`echo ${temp2//,}`
fi

if [ $? == 0 ] 
	then
	echo "step1.1, trimming process sucessful!" >> pipe_processing.log
else 
	echo "step1.1, trimming process fail......" >> pipe_processing.log
	exit 1
fi

# 1.2 fastqc
echo 'fastqc is processing fastq file......'
fastqc -t $threads 'Trimed_'$name*'.fastq' -o . 

if [ $? == 0 ] 
	then
	echo "step1.2, fastqc process sucessful!" >> pipe_processing.log
else 
	echo "step1.2, fastqc process fail......" >> pipe_processing.log
	exit 1
fi

for zip in `ls | grep fastqc.zip`
do	
unzip $zip
rm $zip
mv ${zip%.zip} 'step1.2_'${zip%.zip}
done


# 1.3 fastqc data collection (rely on the output data structure of Fastqc, double-check if it's updated)
echo -e "filename\tdeduplication_percentage\tmarker" > 'dedup_percentage_'$name'.result'
for file in `ls -d *fastqc/`
do
	cd $file
	temp=`echo ${file##Trimed_}`
	out_name=`echo ${temp%*_fastqc/}`
	out_value=`grep 'Total Deduplicated Percentage' fastqc_data.txt | awk '{print $4}'`
	echo -e "$out_name\t$out_value\t$marker" >> ../'dedup_percentage_'$name'.result'
	echo -e "item\t$out_name\t$out_name" > 'duplication_summary_'$out_name'.result'
	grep 'Sequence Duplication Levels' -A 15 fastqc_data.txt >> 'duplication_summary_'$out_name'.result'
	mv 'duplication_summary_'$out_name'.result' ../'data_collection_'$name
	echo -e "$out_name\tfastqc_test" > 'fastqc_summary_'$out_name'.result'
	awk -F "\t" '{print $1,$2}' OFS='\t' summary.txt >> 'fastqc_summary_'$out_name'.result'
	mv 'fastqc_summary_'$out_name'.result' ../'data_collection_'$name
	cd ..
done

if [ $? == 0 ] 
	then
	echo "step1.3, fastqc data_collection process sucessful!" >> pipe_processing.log
else 
	echo "step1.3, fastqc data_collection process fail......" >> pipe_processing.log
fi 

mv 'dedup_percentage_'$name'.result'  ./'data_collection_'$name


# 1.4, get PE data R1 R2 deduplication difference percentage 
if [[ $types == PE ]];
then
per1=`tail -n 2 ./'data_collection_'$name/'dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '1p'`
per2=`tail -n 2 ./'data_collection_'$name/'dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '2p'`
dif=`echo "scale=2; ($per1-$per2)*200/($per1+$per2)" | bc -l`
else
dif=0
fi

if [ $? == 0 ] 
	then
	echo "step1.4, calculate replicate difference process sucessful!" >> pipe_processing.log
else 
	echo "step1.4, calculate replicate difference process fail......" >> pipe_processing.log
fi


###################################################################################################

# step2, alignment and data processing
# 2.1 alignment by HISAT2
# STAR  --runThreadN 24 --genomeDir $star_ref --readFilesIn 'Trimed_'$name*'.fastq' --outSAMtype BAM SortedByCoordinate \
# --outFileNamePrefix 'step2.1_Star_'$name'_'
if [[ $types == PE ]];
then
hisat2 -q -x /home/shaopengliu/resources/mm10/hisat2_index/ht2_idx_mm10 -1 'Trimed_'$name'_1.fastq'  -2 'Trimed_'$name'_2.fastq'  2> step2.1_hisat2_summary.txt | samtools view -bS - | samtools sort - -o step2.1_hisat2_sorted_$name'.bam'
echo "step2.1, mapping as PE" >> pipe_processing.log
else
hisat2 -q -x /home/shaopengliu/resources/mm10/hisat2_index/ht2_idx_mm10 -U 'Trimed_'$name'.fastq' 2> step2.1_hisat2_summary.txt | samtools view -bS - | samtools sort - -o step2.1_hisat2_sorted_$name'.bam'
echo "step2.1, mapping as SE" >> pipe_processing.log
fi

if [ $? == 0 ] 
	then
	echo "step2.1, HISAT2 alignment process sucessful!" >> pipe_processing.log
else 
	echo "step2.1, HISAT2 alignment process fail......" >> pipe_processing.log
	exit 1
fi

# step2.1_hisat2_summary.txt
align_0_time_num=`grep "0 times" step2.1_hisat2_summary.txt | head -1 | awk '{print $1}'`
exact_1_time_num=`grep "exactly 1 time" step2.1_hisat2_summary.txt | head -1 | awk '{print $1}'`
gt_1_time_num=`grep ">1 times" step2.1_hisat2_summary.txt | head -1 | awk '{print $1}'`



for file in `ls *fastq 2> /dev/null`
do
if [[ $file != $R1 ]] && [[ $file != $R2 ]]
then
rm $file 2> /dev/null
fi
done

#index for following step
samtools index step2.1_hisat2_sorted_$name'.bam'

# 3,  QC count: RSeQC
# 3.1, gene coverage:
python2.7 $RSeQC_script'/geneBody_coverage.py' -i step2.1_hisat2_sorted_$name'.bam' -r $ref_gene_bed -o 'step3.1_'$name 
auc=`head -1 step3.1_*.r | awk -F "[(]" '{print $2}' | sed 's/)//' | awk -F "[,]" '{for(i=1;i<=NF;i++) t+=$i; print t}'`

# 3.2, RNA fragment size
python2.7 $RSeQC_script'/RNA_fragment_size.py' -r $ref_gene_bed -i step2.1_hisat2_sorted_$name'.bam'  > 'step3.2_'$name'.fragSize'

# 3.3, inner distance
python2.7 $RSeQC_script'/inner_distance.py' -i step2.1_hisat2_sorted_$name'.bam' -o 'step3.3_'$name -r $ref_gene_bed

# 3.4, insertion profile
python2.7 $RSeQC_script'/insertion_profile.py' -i step2.1_hisat2_sorted_$name'.bam' -s $types -o 'step3.4_'$name

# 3.5, read GC
python2.7 $RSeQC_script'/read_GC.py' -i step2.1_hisat2_sorted_$name'.bam' -o 'step3.5_'$name

# 3.6, reads distribution
python2.7 $RSeQC_script'/read_distribution.py' -i step2.1_hisat2_sorted_$name'.bam'  -r $ref_gene_bed > 'step3.6_'$name'_read_distribution.txt'

#3.7, reads duplication
python2.7 $RSeQC_script'/read_duplication.py' -i step2.1_hisat2_sorted_$name'.bam' -o 'step3.7_'$name

#3.8, junction annotation
python2.7 $RSeQC_script'/junction_annotation.py' -i step2.1_hisat2_sorted_$name'.bam'  -o 'step3.8_'$name -r $ref_gene_bed

#3.9, RPKM_saturation.py
python2.7 $RSeQC_script'/RPKM_saturation.py' -i step2.1_hisat2_sorted_$name'.bam'  -o 'step3.9_'$name -r $ref_gene_bed

# 3.10, preseq
# bamToBed -i 'Star_'$name'_Aligned.sortedByCoord.out.bam' > temp.bed
# $preseq lc_extrap -o 'yield_'$name'.result' temp.bed  
$preseq lc_extrap -o 'step3.10_yield_'$name'.result'  -B step2.1_hisat2_sorted_$name'.bam'
cp 'step3.10_yield_'$name'.result'  ./data_collection_$name/'yield_'$name'.result'



# 4.1, feature count
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
	featureCounts -a $annotation_file -T $threads -O  \
	 -o step4.1_gene_name_fc_$name  -g gene_name step2.1_hisat2_sorted_$name'.bam'
	featureCounts -a $annotation_file  -T $threads -O  \
	 -o step4.1_gene_type_fc_$name -g gene_type step2.1_hisat2_sorted_$name'.bam'
fi

# 4.2, FC collection
echo -e "gene_name\tlength\tfragment_count" > 'step4.2_gene_name_count_'$name'.txt'
tail -n +3 step4.1_gene_name_fc_$name | awk '{print $1"\t"$6"\t"$7}' | sort -k3,3rn >> 'step4.2_gene_name_count_'$name'.txt'

echo -e "gene_type\tlength\tfragment_count" > 'step4.2_gene_type_count_'$name'.txt'
tail -n +3 step4.1_gene_type_fc_$name | awk '{print $1"\t"$6"\t"$7}' | sort -k3,3rn >> 'step4.2_gene_type_count_'$name'.txt'

cp step4.2*txt  ./data_collection_$name
awk '{print $0}' step3.5_*.GC.xls  > './data_collection_'$name'/step3.5_'$name'.GC.txt'
awk '$0~/====/ {p=1}; p;' step3.6_*read_distribution.txt | tail -n +2 | head -n -1 | awk '{print $1"\t"$2"\t"$3"\t"$4}' > step3.6_$name'_reads_distri_table.txt'
mv step3.6_$name'_reads_distri_table.txt'   ./data_collection_$name

# ratio of reads with feature
total_cc=`awk '{s+=$2}END{print s}' step4.1_gene_name_fc*summary`
reads_with_feature=`grep Assigned step4.1_gene_name_fc*summary | awk '{print $2}'`
rates_with_feature=`python -c "print(1.0*$reads_with_feature/$total_cc)"`

cp step4.2_gene_name*txt ./data_collection_$name
total_count=`awk '{s+=$3}END{print s}' step4.2_gene_name*txt`
thres=`python -c "print(1.0* $total_count/1000000)"`
exp_gene_num=`awk -v thres=$thres '$3>thres' step4.2_gene_name*txt | wc -l`
#summarize together:
echo -e "name\traw_reads\tremoved_reads\talign_0_time\texactly_1_time\tgt_1_time\treads_ratio_with_gene_feature\tauc_genebody\texp_gene_num" > RNA-seq_qc_collection.txt
echo -e "$name\t$raw_reads\t$removed_reads\t$align_0_time_num\t$exact_1_time_num\t$gt_1_time_num\t$rates_with_feature\t$auc\t$exp_gene_num" >> RNA-seq_qc_collection.txt
mv RNA-seq_qc_collection.txt ./data_collection_$name

# 5, multiqc, v1.3+ (need)
multiqc .

# clean results
rm refined_chrom_size.txt 

# 6, get summarize table and json output
cd  ./data_collection_$name
Rscript /home/shaopengliu/test_rna/test_pipe/summarize.R  $name $species
sed 's/\[/{/g' $name'_report.json' | sed '/    {/d' | sed '/\]/d' | sed 's/    }/  },/g' | tac | sed '2s/},/}/g' | tac > $name'.json'
rm $name'_report.json'
cd ..

echo "processing finished"
date