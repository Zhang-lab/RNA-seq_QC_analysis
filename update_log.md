### Current version: `190321`       
Last update: 03/21/2019  
  
### Update record:  
03/21/2019 TaRGET & HPC  
1. fix some typo on indexing step of STAR bedgraph file  
2. fix mislabel on gtf file annotation in json file  

03/11/2019 TaRGET & HPC  
1. revise the logic for strand determination; previously it will choose the strand that yield highest output, but recently our finding on TaRGET stranded data show many cases where the s0 is close to s1 / s2, which may yield wrong estimate. So current logic is that given a stranded data, it shall has signalicant higher assigned reads than its counterpart (we use 2-folder).  

02/27/2019 TaRGET & HPC  
1. replace mm10 reference from gtf v15 -> v20  

12/07/2018 TaRGET  
1. adjust the threshold for passing "uniquely mapped ratio" from 0.65 to 0.35  

10/12/2018 TaRGET  
1. add target specific version target_181012:  
    1) add score matrix 1 / 0 for 3 QC term  
    2) add preseq file check for json file  
2. add target specific integrative analysis version target_integrative:  
    1) hard trim all data into 50bp from 3' end  

09/28/2018 v4  
1. add featureCount automatic strand check for input data, now "-s" is NOT an optional input parameter anymore  
2. explicitely specify -F "\t" for awk, otherwise it may misunderstand danRer10 fc output  
3. details on STAR parameter choice:  
    1) generating signal with --outWigType requires sorted BAM (by coordinate)  
    2) featureCounts requires sort by name (default STAR output)  


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
