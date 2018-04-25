date

for file in `find . -name "Processed*"`
do
name=`echo ${file#*Processed_}`
cd $file
echo "processing $name"

total_cc=`grep "Uniquely mapped reads number" step2.1_Star_$name'_Log.final.out' | awk '{print $6}'`
reads_with_feature=`grep Assigned step4.1_gene_name_fc*summary | awk '{print $2}'`
rates_with_feature=`python -c "print(1.0*$reads_with_feature/$total_cc)"`

cat ./data_collection_$name'/RNA-seq_qc_collection.txt' | awk -v rates=$rates_with_feature 'NR==2{$7=rates};1' OFS="\t"  > temp.txt \
    && mv temp.txt  ./data_collection_$name'/RNA-seq_qc_collection.txt'

cd -
done

echo "whole pipe done"
date






