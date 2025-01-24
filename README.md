# TriTry-ncRNA
TriTry-ncRNA is a pipeline in order to identify putative non-Coding RNA in Trypanosomatids organism based on manages transcripts from RNA-seq-data. This tool begins with the detection of the set of transcripts mapped by sequencing the total RNA of the three life stages of Leishmania braziliensis: procyclic promastigote (PRO), metacyclic promastigote (META) and amastigote (AMA).

## Introduction
* 2_identify_transcript.py: Identify in all chromossomal positions the coverage equal & up of 50x or 100x reads for strand + & - strand.
* 3_identify_possible_ncRNA_lncRNA.py: Identify lnc-RNA (>200pb) with >= 50x cov and snc-RNA (<200 pb) with >= 100x cov.
* 4_annotation_ncRNA_lncRNA.py: Identify the location of non-coding RNA in coding or intergenic regions. 
* 5_identify_overlap_nc-lncRNA.py: Identify the snc-RNA overlapping lnc-RNA.
* 6_identify_ptu.py: To identify: PTU regions.
* 7_parser_UTR_sense_antisense.py: Add information of sense of non-coding RNA in relation from PTUs and UTR regions.
* 8_filterpfam.py: Filter non-coding RNA from blast against pfam proteins.
* 9_select_ncrna.py: Final selection by directory.
* 10_remake_output.py: Create gff and bed format output files.

The pipeline process all sub-modules by directory, each directory 
can reprsent a life cycle, that contain the biological repetions in fastq format (*_1/2.fastq.gz). Posteriolly this are merge for all total non-coding RNAs.

##Input files
TriTry-ncRNA users should have the following minimal input files:
- lisdir.txt (-dir_list) : File contain a list of
 directories path that contain all fastq files.
- genome.fasta (-fasta): Reference genome FASTA file.
- genome.gff (-gff): Annotation GFF file.

## Instalation and requirements
The following software and libraries must be installed on your machine:
* Fastqc
* Bowtie2 
* Samtools
* Picard
* Bedtools
* BLAST+ 
* Python3

## Installation
Clone the TriTry-ncRNA git:
```
git clone https://github.com/AKCLAB/TriTry-ncRNA.git
```

Go to scripts directory
```
cd TriTry-ncRNA
```

## Test with example:
```
cd test/
```




