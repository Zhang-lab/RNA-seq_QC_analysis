#!/bin/bash
# pipe usage:
# user@domain: path_to_pipe/pipe.sh -g <hg38/mm10/danRer10>  -r <PE/SE>  -o file1  -p file2

# pipe start
###################################################################################################
# Preparation:
date
pipe_version='targetv2'

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
md5=`md5sum $0 | awk '{print $1}'`

# read parameters
while getopts m:t:g:o:p:r:a:b:f  opts
do case "$opts" in
m) marker="$OPTARG";;	# default 'unmarked'
t) threads="$OPTARG";;	# default 24
g) species="$OPTARG";;	# hg19, hg38, mm9, mm10, danRer10
o) R1="$OPTARG";;    # PE read 1, or the SE file, or the sra file
p) R2="$OPTARG";;    # PE read 2. 
a) adapter_1="$OPTARG";;   # add adapter1
b) adapter_2="$OPTARG";;    # add adapter2
r) types="$OPTARG";;	# PE or SE;
f) fast_mode=1;;	# fast mode
h) echo "usage:  path-to-pipe/pipe.sh  -g <hg38/hg19/mm10/mm9/danRer10/personal>  -r <PE/SE> -o read_file1  -p read_file2 (if necessary)"
exit;;
[?]) "Usage: ./pipe.sh  -g <hg38/hg19/mm10/mm9/danRer10> -o <read_file1>  -r <PE/SE>";;
esac
done

if [ -z "$fast_mode" ]
then
echo "  "
echo "processing whole pipe"
else
	echo "  "
	echo "choosing fast mode, skipping multiQC steps!!!"
fi


if [ -z "$threads" ]
then
threads=24
fi

if [ -z "$marker" ]
then
echo "you didn't specify the marker, are you sure to keep the default unmarked"
marker='unmarked'
fi

if [ -z "$adapter_1" ]
	then
	adapter_1="AGATCGGAAGAGCACACGTCTGAAC"
	adapter_2="AGATCGGAAGAGCGTCGTGTAGGGA"
	notrim="--no-trim"
	echo "No adapter trimming happened because there is no adapter input, only quality trimming on reads"
fi


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
ln -rs $R1    ./'Processed_'$name/
ln -rs $raw1  ./'Processed_'$name/  2> /dev/null
ln -rs $raw2  ./'Processed_'$name/  2> /dev/null
cd ./'Processed_'$name/
mkdir 'data_collection_'$name
touch pipe_processing.log
source $pipe_path'/rna_seq_qc_source.sh' $species

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
	cutadapt $notrim -j $threads -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=36  -o 'Trimed_'$name'_1.fastq' -p  'Trimed_'$name'_2.fastq'  $raw1 $raw2  > 'step1.1_'$name'_cutadapt_PE.trimlog'  
	temp=`grep "Total read pairs processed:" step1.1_*trimlog | awk '{print $5}'`  
	raw_reads=`echo ${temp//,}`  
	temp2=`grep "Pairs written" step1.1_*trimlog | awk '{print $5}'`
	written_reads=`echo ${temp2//,}`
elif [[ $types == SE ]];
	then
	echo 'trimming RNA-seq SE reads by cutadapt'
	cutadapt $notrim -j $threads -a $adapter_1 --quality-cutoff=15,10 --minimum-length=36  -o  'Trimed_'$name'.fastq' $raw1 > 'step1.1_'$name'_cutadapt_SE.trimlog'  
	temp=`grep "Total reads processed:" step1.1_*trimlog | awk '{print $4}'`
	raw_reads=`echo ${temp//,}` 
	temp2=`grep "Reads written" step1.1_*trimlog | awk '{print $5}'`
	written_reads=`echo ${temp2//,}`
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
# 2.1 alignment by STAR v2.5.4b
STAR  --runThreadN $threads --genomeDir $star_ref --readFilesIn 'Trimed_'$name*'.fastq'  --outFileNamePrefix  'step2.1_Star_'$name'_'  --outWigType bedGraph --outWigNorm RPM --outSAMtype BAM  SortedByCoordinate  


if [ $? == 0 ] 
	then
	echo "step2.1, STAR alignment process sucessful!" >> pipe_processing.log
else 
	echo "step2.1, STAR alignment process fail......" >> pipe_processing.log
	exit 1
fi

awk '{print $1"\t"$2"\t"$3"\t-"$4}' 'step2.1_Star_'$name'_Signal.Unique.str2.out.bg' > test.symbol \
    && cat  test.symbol  'step2.1_Star_'$name'_Signal.Unique.str1.out.bg' > 'step2.1_Star_'$name'.Unique.stranded.bg' \
    && rm test.symbol \
    && bedSort 'step2.1_Star_'$name'.Unique.stranded.bg' 'step2.1_Star_'$name'.Unique.stranded.bg' \
    && bgzip 'step2.1_Star_'$name'.Unique.stranded.bg' \
    && tabix -p bed 'step2.1_Star_'$name'.Unique.stranded.bg.gz'

# step2.1_STAR_summary.txt
input_reads=`grep "Number of input reads" step2.1_Star_$name'_Log.final.out' | awk -F "|" '{print $2}'`
exact_1_time_num=`grep "Uniquely mapped reads number"  step2.1_Star_$name'_Log.final.out' | awk '{print $6}'`
gt_1_time_num=`grep "Number of reads mapped to" step2.1_Star_$name'_Log.final.out'  | awk -F "|" '{print $2}' | awk '{s+=$1}END{print s}'`
align_0_time_num=`python -c "print($input_reads-$exact_1_time_num-$gt_1_time_num)"`

bam_file='step2.1_Star_'$name'_Aligned.sortedByCoord.out.bam'

cal_qc () {
#index for following step
samtools index $bam_file

# 3.10, preseq
awk '{ if (NR%4==2) print substr($0,1,20); }'  `ls 'Trimed_'$name*'.fastq' | head -1`  | sort | uniq -c | awk '{ print $1 }' > counts.txt  && \
$preseq lc_extrap -x 100 -v -Q -o 'step3.10_yield_'$name'.result'  -V counts.txt && \
cp 'step3.10_yield_'$name'.result'  ./data_collection_$name/'yield_'$name'.result'  && \
rm counts.txt 

# 3.1, gene coverage:
$RSeQC_python $RSeQC_script'/geneBody_coverage.py' -i $bam_file -r $sub_ref_gene -o 'step3.1_'$name  && \
auc=`head -1 step3.1_$name'.geneBodyCoverage.r' | awk -F "[(]" '{print $2}' | sed 's/)//' | awk -F "[,]" '{for(i=1;i<=NF;i++) t+=$i; print t}'`  && \
cp step3.1_$name'.geneBodyCoverage.txt'  ./data_collection_*  

# 3.2, RNA fragment size
$RSeQC_python $RSeQC_script'/RNA_fragment_size.py' -r $sub_ref_gene -i $bam_file  > 'step3.2_'$name'.fragSize' 

# 3.3, inner distance
$RSeQC_python $RSeQC_script'/inner_distance.py' -i $bam_file -o 'step3.3_'$name -r $sub_ref_gene  

# 3.5, read GC
$RSeQC_python $RSeQC_script'/read_GC.py' -i $bam_file -o 'step3.5_'$name  

# 3.6, reads distribution
$RSeQC_python $RSeQC_script'/read_distribution.py' -i $bam_file  -r $sub_ref_gene > 'step3.6_'$name'_read_distribution.txt'  

#3.8, junction annotation
#keep full gene list for this one:
$RSeQC_python $RSeQC_script'/junction_annotation.py' -i $bam_file  -o 'step3.8_'$name -r $ref_gene_bed  && \
sed 1d step3.8_$name'.junction.xls' | cut -f 5 | sort | uniq -c > a.txt  && \
sed 1d step3.8_$name'.junction.xls' | cut -f 4,5 | awk '{if ($2=="complete_novel") s2+=$1; else if ($2=="annotated") s1+=$1; else if ($2=="partial_novel") s3+=$1}END{print s1"\n"s2"\n"s3}' | paste <(awk '{print $2,$1}' OFS="\t" a.txt) -  > step3.8_$name'.result'  && \
mv step3.8_$name'.result' ./data_collection*  
}


if [ -z "$fast_mode" ]
then
cal_qc
fi

# 4.1, feature count
if [[ $types == PE ]];
	then
	echo 'featureCounts on PE data'
	featureCounts -a $annotation_file -p -T 4 -O  \
	 -o step4.1_gene_name_fc_$name  -g gene_name $bam_file  
	featureCounts -a $annotation_file -p -T 4 -O  \
	 -o step4.1_gene_type_fc_$name -g gene_type $bam_file   
	featureCounts -a $annotation_file -p -T 4 -O  \
	 -o step4.1_transcript_fc_$name -g transcript_id $bam_file   
elif [[ $types == SE ]];
	then
	echo 'featureCounts on SE data'
	featureCounts -a $annotation_file -T 4 -O  \
	 -o step4.1_gene_name_fc_$name  -g gene_name $bam_file  
	featureCounts -a $annotation_file  -T 4 -O  \
	 -o step4.1_gene_type_fc_$name -g gene_type $bam_file  
	featureCounts -a $annotation_file  -T 4 -O  \
	 -o step4.1_transcript_fc_$name -g transcript_id $bam_file  
fi



# 4.2, FC collection
echo -e "gene_name\tlength\tfragment_count" > 'step4.2_gene_name_count_'$name'.txt'
tail -n +3 step4.1_gene_name_fc_$name | awk '{print $1"\t"$6"\t"$7}' | sort -k3,3rn >> 'step4.2_gene_name_count_'$name'.txt'

echo -e "gene_type\tlength\tfragment_count" > 'step4.2_gene_type_count_'$name'.txt'
tail -n +3 step4.1_gene_type_fc_$name | awk '{print $1"\t"$6"\t"$7}' | sort -k3,3rn >> 'step4.2_gene_type_count_'$name'.txt'

echo -e "transcript_id\tlength\tfragment_count" > 'step4.2_transcript_count_'$name'.txt'
tail -n +3 step4.1_transcript_fc_$name | awk '{print $1"\t"$6"\t"$7}' | sort -k3,3rn >>  'step4.2_transcript_count_'$name'.txt'

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
uniq_ratio=`python -c "print(1.0* $exact_1_time_num/$written_reads)"`
#summarize together:
echo -e "name\traw_reads\twritten_reads\talign_0_time\texactly_1_time\tuniq_ratio\tgt_1_time\treads_ratio_with_gene_feature\tauc_genebody\texp_gene_num" > RNA-seq_qc_collection.txt
echo -e "$name\t$raw_reads\t$written_reads\t$align_0_time_num\t$exact_1_time_num\t$uniq_ratio\t$gt_1_time_num\t$rates_with_feature\t$auc\t$exp_gene_num" >> RNA-seq_qc_collection.txt
mv RNA-seq_qc_collection.txt ./data_collection_$name

# 4.3, FC saturation analysis
cut -f 1,3 'step4.2_gene_name_count_'$name'.txt' | sed 1d | awk '$2>0' > temp.txt
while read i  # 3min for a huge file
do
	yes `echo $i | cut  -d " " -f 1` |  head -`echo $i | cut  -d " " -f 2` >> fc_subsample.txt
done < temp.txt

all_fc=`wc -l fc_subsample.txt | awk '{print $1}'`
awk -v sum_count=$all_fc '{print $1,$2*1000000/sum_count}' OFS="\t" temp.txt | sort > raw_CPM.txt && \
	rm temp.txt

cat fc_subsample.txt |  sort | uniq -c | awk -v sample_ratio=$all_fc '{print $2,$1*1000000/sample_ratio}' OFS="\t" | awk '$2>=1' >  raw_CPM.txt
lt10=`awk '$2<10' raw_CPM.txt | wc -l`
bt50=`awk '$2>=10 && $2<50' raw_CPM.txt | wc -l`
gt50=`awk '$2>=50' raw_CPM.txt | wc -l`
all_gene=`cat raw_CPM.txt | wc -l`

# subsampling:
for number in 10 20 30 40 50 60 70 80 90
do
	sample_ratio=$((all_fc * $number / 100 ))
	echo "$sample_ratio out of total $all_fc"
	shuf fc_subsample.txt | head -$sample_ratio | sort | uniq -c | awk -v sample_ratio=$sample_ratio '{print $2,$1*1000000/sample_ratio}' OFS="\t" > b.txt  && \
	join -a 1 -j 1 raw_CPM.txt b.txt | awk 'sqrt(($2-$3)^2)/$2>0.1' > c.txt
	echo -e "`awk '$2<10' c.txt | wc -l`\t`awk '$2>=10 && $2<50' c.txt | wc -l`\t`awk '$2>=50' c.txt | wc -l`\t`cat c.txt | wc -l`" >> step4.3_CPM_saturation_$name'.txt'
done
awk -v lt10=$lt10 -v bt50=$bt50 -v gt50=$gt50 -v all_gene=$all_gene '{print lt10-$1, bt50-$2, gt50-$3, all_gene-$4 }' OFS="\t"  step4.3_CPM_saturation_$name'.txt' > temp.txt && mv temp.txt   step4.3_CPM_saturation_$name'.txt'
echo -e "$lt10\t$bt50\t$gt50\t$all_gene" >> step4.3_CPM_saturation_$name'.txt' && \
	paste <(seq 10 10 100) step4.3_CPM_saturation_$name'.txt' | cat <(echo -e "percent\tlt10\tbt50\tgt50\ttotal") -> temp.txt && mv temp.txt  step4.3_CPM_saturation_$name'.txt' 

rm a.txt b.txt c.txt
cp step4.3_CPM_saturation_$name'.txt' ./data_collection_*


# 5, multiqc, v1.3+ 
multiqc .

# clean results
rm refined_chrom_size.txt 
rm fc_subsample.txt 
rm log.txt 

for file in `ls *fastq 2> /dev/null`
do
if [[ $file != $R1 ]] && [[ $file != $R2 ]]
then
rm $file 2> /dev/null
fi
done

# 6, get summarize table and json output
time=`head -1 pipe_processing.log | sed 's/ /_/g'`
host="zhanglab/rna-seq targetv2"
image_id=`bash $pipe_path'/find_image_ID_digest.sh' $host  2> /dev/null | awk '{print $2}'`
if [ -z "$image_id" ]
    then
    image_id="failed_to_get_id"
fi

cd  ./data_collection_$name
Rscript  $pipe_path'/summarize.R'  $name $species  $pipe_version $time $image_id
sed 's/\[/{/g' $name'_report.json' | sed '/    {/d' | sed '/\]/d' | sed 's/    }/  },/g' | sed 's/"!/{/g' | sed 's/!"/}/g' | sed 's/"?/[/g' | sed 's/?"/]/g' | sed 's/@/"/g' | tac | sed '2s/},/}/g' | tac | sed "s/MD5ToBeChange/$md5/g" > $name'.json'
rm $name'_report.json'
cd ..

echo "processing finished"
date

