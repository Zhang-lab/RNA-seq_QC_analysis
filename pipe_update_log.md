### Current version: v1.1b  
Last update: 03/28/2018  
  
### Update record: 
03/28/2018 v1.1b  
1. build STAR ref by version 2.5.4b, gencode annotation vM15  
2. change STAR alignment code (stranded)  
3. reduce gene list for RSeQC RPKM saturation calculation, randomly sampling 5k genes from refseq gene list)  
4. add cutadapt filter: if no adapter input there would be only quality trim  


03/16/2018 v1.1a  
1. add --outBAMsortingThreadN parameter into STAR  
2. adjust STAR threads setting from 24 to $threads  
