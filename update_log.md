### Current version: `v4`      
Last update: 09/05/2018  
  
### Update record:  
09/05/2018 v4  
1. remove strand option for featureCount (added in v2.2); now we will use chr10 reads to do a test to determine the strand information  
2. for mm10, add globin gene counts from featureCount into json output  

07/13/2018 v3.1  
1. adjust json output  

06/06/2018 v3  
1. add RSEM for transcript level quantification (featureCounts not as expected)  
2. modify output json files  
3. modify output folder structure and file names  
4. verify strand information on STAR alignment step  

05/15/2018 v2.2  
1. remove -O option in featureCount (assign multiple features)  
2. add -Q 10 for featureCount  
3. add strand option for featureCount, now -s is a required option  

05/10/2018 v2.1  
1. add merge bedgraph step for STAR output  

05/07/2018 targetv2  
1. fix version for target (v2)  
2. fix docker image and image id for now  

04/12/2018 v1.2  
1. update cutadapt from v1.14 -> v1.16, now support multi-core trimming  


03/28/2018 v1.1b  
1. build STAR ref by version 2.5.4b, gencode annotation vM15  
2. change STAR alignment code (stranded)  
3. reduce gene list for RSeQC RPKM saturation calculation, randomly sampling 5k genes from refseq gene list  
4. add cutadapt filter: if no adapter input there would be only quality trim  


03/16/2018 v1.1a  
1. add --outBAMsortingThreadN parameter into STAR  
2. adjust STAR threads setting from 24 to $threads  
