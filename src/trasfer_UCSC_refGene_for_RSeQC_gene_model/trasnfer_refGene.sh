#!/bin/bash

input_file=$1
echo "processing $input_file"

awk -F "\t" '{print $3,$5,$6,$2,0,$4,$7,$8,$12,$9,$10,$11}' OFS="\t"   $input_file > temp.bed

cat <<EOF > temp_modify.R
read.table("temp.bed", sep="\t", stringsAsFactors = FALSE) -> t1
old_col_2=as.vector(t1\$V2)
old_col_11=strsplit(t1\$V11, ",")
old_col_12=strsplit(t1\$V12, ",")

new_col_11<-vector()
new_col_12<-vector()

for ( i in 1:dim(t1)[1] ) {
    new_col_12[i]=paste( paste((as.numeric(as.character(unlist(old_col_11[i]))) - old_col_2[i]), collapse=","), "," , sep="")
    new_col_11[i]=paste( paste((as.numeric(as.character(unlist(old_col_12[i]))) - as.numeric(as.character(unlist(old_col_11[i]))) + 1), collapse=","), "," , sep="")
}

out<-t1
out\$V11=new_col_11
out\$V12=new_col_12

write.table(out, file="modified_temp.bed", quote=F, sep="\t",col.names=F,row.names=F)
EOF

Rscript temp_modify.R && rm temp_modify.R temp.bed
echo "the output file is modified_temp.bed, please rename it to your own file"
date


