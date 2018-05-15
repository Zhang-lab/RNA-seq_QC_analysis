#!/bin/bash
date_now=`date +"%m-%d-%y"`
mkdir rna_qc_collection_$date_now
cd rna_qc_collection_$date_now
temp=`pwd`

#data collection
for file in `find .. -name "data_collection_*"`   
do
	name=`echo ${file#*data_collection_}`
	#1, gene type distribution
	cd $file
	cp RNA-seq_qc_collection.txt  $temp'/RNA_qc_col_'$name
	cp step4.2_gene_type_count*txt $temp
	cp step4.2_gene_name_count*txt $temp
	# for old pipe:
	#awk '$0~/====/ {p=1}; p;' step3.6_*read_distribution.txt | tail -n +2 | head -n -1 | awk '{print $1"\t"$2"\t"$3"\t"$4}' > step3.6_$name'_reads_distri_table.txt'
	cp step3.6_$name'_reads_distri_table.txt' $temp
	cd $temp
done

#1, merge RNA-qc table:
head -1 `ls RNA_qc_col* | head -1` > merged_RNA_qc_table.txt
for file2 in `ls RNA_qc_col*`
do
	sed -n '2p' $file2 >> merged_RNA_qc_table.txt
	rm $file2
done

#2, merge gene-type count:
tail -n +2 `ls step4.2_gene_type_count_*txt | head -1` | sort | awk '{print $1"\t"$2}' | cat <(echo -e "type\tlength") - > temp.txt
for file3 in `ls step4.2_gene_type_count_*txt`
do
	f3_name=`echo ${file3#step4.2_gene_type_count_}`
	f3_name2=`echo ${f3_name%.txt}`
	tail -n +2 $file3 | sort | awk '{print $3}' | cat <(echo "$f3_name2") -   | paste temp.txt - > temp2.txt
	mv temp2.txt temp.txt
	rm $file3
done
mv temp.txt merged_gene_type_count.txt

#3, merge reads distribution
cat `ls step3.6*reads_distri_table.txt | head -1` | awk '{print $1"\t"$2}' > merged_reads_distri_in_genome.txt
for file4 in `ls step3.6*reads_distri_table.txt`
do
	f4_name=`echo ${file4#step3.6_}`
	f4_name2=`echo ${f4_name%_reads_distri_table.txt}`
	tail -n +2 $file4 | awk '{print $3}' | cat <(echo "tag_count_$f4_name2") - | paste merged_reads_distri_in_genome.txt - > temp.txt
	mv temp.txt merged_reads_distri_in_genome.txt
	rm $file4
done

#4, merge expression table
tail -n +2 `ls step4.2_gene_name_count_*txt | head -1` | sort | awk '{print $1"\t"$2}' |  cat <(echo -e "name\tlength") - > temp.txt
for file4 in `ls step4.2_gene_name_count_*txt`
do
	f4_name=`echo ${file4#step4.2_gene_name_count_}`
	f4_name2=`echo ${f4_name%.txt}`
	tail -n +2 $file4 | sort | awk '{print $3}'  | cat <(echo "$f4_name2") - | paste temp.txt - > temp2.txt
	rm $file4
	mv temp2.txt temp.txt
done
mv temp.txt  merged_gene_name_expression.txt








