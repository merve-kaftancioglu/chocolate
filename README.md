# chocolate
ChIP-seq analysis pipeline (in progress)

# Overview

The **chocolate** pipeline accepts fastq files, and emits various QC files (for reads, peaks, profiles), called peaks, and differentially bound peaks.  

The pipeline is based on an existing pipline by NF-Core (https://github.com/nf-core/chipseq), as well as a pipeline that I created previously. At the moment, the pipeline assumes single-end sequencing, but can be extended to paired-end sequencing later.

**chocolate** performs the following steps:

_Steps_
 - Read quality asessment (`FastQC`, `FastQScreen`)
 - Read trimming (`bbduk`)
 - Read alignment (`BWA mem`)
 - Alignment file sorting and indexing (`samtools sort` and `index`)
 - Marking duplicate reads (`Picard MarkDuplicates`)
 - Filter reads (`samtools view`)
 - `-F 4` - exclude non-mapping reads
 - `-F 256` - exclude multi-mapping reads
 - `-F 1024` - exclude PCR or optical duplicates
 - Sort and index (`samtools`)
 - ChIP QC steps for enrichment (`deepTools`, `phantompeakqualtools`, `calculateNSCRSC.r`, `ngsplot`)
 - Call peaks (`MACS2`)
 - ChIP peak annotation (`chippeakanno`)

**chocolate** is implemented in Snakemake and makes use of several bioinformatic tools. Two main files (`chocolate.smk` and `config.yaml`) and a host of supporting files/directories (`rules/`, `snakemake_env.yml`) are necessary. See below. 
