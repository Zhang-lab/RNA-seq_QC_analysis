# RNA_seq_analysis_pipe  
This is for the QC matrix construction, data analysis and visualization for RNA-seq data.  
Current version: `v4`   

Advisor: Bo Zhang  
Contributor: Cheng Lyu and Shaopeng Liu  

For any question, please contact Wustl.Zhanglab@gmail.com  


## Usage: 
Singularity 2-step solution (easiest way)  

Step1. download singularity container (you only need download the containcer for once, then you can use them directly):  
####  
```bash
# download image from local server:  
wget http://regmedsrv1.wustl.edu/Public_SPACE/shaopengliu/Singularity_image/rna-seq/rna-seq_mm10_v4.simg  
```

or for human

```bash
# download image from local server:
wget http://regmedsrv1.wustl.edu/Public_SPACE/shaopengliu/Singularity_image/rna-seq/hg38_rna-seq.simg
```

Step2. process data by the singularity image: 
#### Please run at same directory with your data OR the soft link of your data    
```bash
singularity run -H ./:/scratch rna-seq_mm10_v4.simg -r <SE/PE> -g <mm10>  -o <read_file1>  -p <read_file2>    
```

That's it!

#parameters:  
`-h`: help information  
`-r`: SE for single-end, PE for paired-end  
`-g`: genome reference, one simg is designed for ONLY one species due to the file size. For now the supported genoms are: <mm10/mm9/hg19/hg38/danRer10> (only mm10 in the example).  
`-o`: reads file 1 or the SE reads file, must be ended by .fastq or .fastq.gz or .sra (for both SE and PE)  
`-p`: reads file 2 if input PE data, must be ended by .fastq or .fastq.gz  
`-a`: ADAPT1 for cutadapt  
`-b`: ADAPT2 for cutadapt, if not specified there will be ONLY quality trimming rather than adapter trimming    
`-t`: (optional) specify number of threads to use, default 24  

e.g:
a) mm10 SE data A.fastq  
```bash
singularity run -H ./:/scratch <path_to_simg>  -r SE -g mm10 -o A.fastq  
```
b) hg38 PE data B_1.fastq B_2.fastq  
```bash
singularity run -H ./:/scratch <path_to_simg>  -r PE -g hg38 -o B_1.fastq  -p B_2.fastq  
```
c) danRer10 PE data in sra file C.sra  
```bash
singularity run -H ./:/scratch <path_to_simg> -r PE -g danRer10 -o C.sra  
```


