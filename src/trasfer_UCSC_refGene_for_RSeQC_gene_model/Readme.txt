# last update: 20180927

# RSeQC requires gene model for gene annotation, which can be downloaded from UCSC
# RSeQC only provide hg38, mm10 and DanRer7
# this code is to tranfer downloaded refGene to required format:

# step1, reorganize columns:
UCSC format:
table genePredExt
"A gene prediction with some additional info."
    (
    string name;            "Name of gene (usually transcript_id from GTF)"
    string chrom;           "Chromosome name"
    char[1] strand;         "+ or - for strand"
    uint txStart;           "Transcription start position"
    uint txEnd;             "Transcription end position"
    uint cdsStart;          "Coding region start"
    uint cdsEnd;            "Coding region end"
    uint exonCount;         "Number of exons"
    uint[exonCount] exonStarts; "Exon start positions"
    uint[exonCount] exonEnds;   "Exon end positions"
    uint id;                "Unique identifier"
    string name2;           "Alternate name (e.g. gene_id from GTF)"
    string cdsStartStat;     "enum('none','unk','incmpl','cmpl')"
    string cdsEndStat;       "enum('none','unk','incmpl','cmpl')"
    lstring exonFrames;     "Exon frame offsets {0,1,2}"
    )
* there is one extra column of number at the very beginning



RSeQC required: 
# here are 2 unit_id column
# in the RSeQC provided files, all those 2 columns are all 0, I assume those 2 are designed for different purpose (gene / transcrip?) or just for annotation purpose; simply keep them here

chr	start	end	name	unit_id_or_0	strand	cds_start	cds_end	unit_id	exon_count	exon_length	exon_start_position_to_start_point




