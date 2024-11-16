## CAB small RNA-seq Nextflow pipeline [V1.2.0] based on sRNAtools

## CAB internal wikipage
https://wiki.stjude.org/display/CABI/CAB+miRNAseq+nextflow+pipeline

## Dependencies 
* databases: /research/groups/cab/projects/automapper/common/tchang1/Softwares/sRNAtools/db/  (defined by DBCONFIG.txt) 
* bowtie-1.2.2 index: /research/groups/cab/projects/automapper/common/hjin/bin/bowtie-1.2.2/indexes
* nextflow/21.10.5
* modules: R/4.3.2,perl/5.20.1, samtools/1.9 , bedtools/2.25.0, bowtie/1.2.2, trimgalore/0.6.6,fastqc/0.11.8, python37, 

## Input

* format1 (2 columns) for single_end_oneLane: 
  SampleName R1.fastq.gz 
* format2 (2 columns) for single_end_multiLanes: 
  SampleName L1_R1.fastq.gz,L2_R1.fastq.gz 

## Run - example1
```bash
export LD_LIBRARY_PATH="/hpcf/authorized_apps/rhel8_apps/gcc/13.1.0/install/lib64:$LD_LIBRARY_PATH"
module nextflow/21.10.5
nextflow run sRNAtools.nf --help 


Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/research/rgs01/scratch_lsf/java -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1
N E X T F L O W  ~  version 21.10.6
Launching `./sRNAtools.nf` [silly_spence] - revision: 7d0f8d63c3

        Welocme to run Nextflow Pipeline sRNAtools.nf [version 1.3.0, 11/15/2024]
       Usage:
        A typical command for running the pipeline is as follows:
        nextflow run sRNAtools.nf -profile local --fqlist fq1.lst --Trim_Mode 1 --outdir run1 --prefix hendegrpq
        nextflow run sRNAtools.nf -profile local --fqlist fq2.lst --Trim_Mode 2 --outdir run2 --prefix sapkogrp

        Configurable arguments:
        --outdir                      The directory for the pipeline output
        --prefix                      output prefix for your analysis
        --A                           adapter sequence for read 1 if Trim_Mode is 3
        --fqlist                      Fastq file list that you want to call sRNAtools
                 format1 (2 columns) for single_end_oneLane:
                          SampleName L001_R1_001.fastq.gz
                 format2 (2 columns) for single_end_multiLanes:
                          SampleName L001_R1_001.fastq.gz,L002_R1_001.fastq.gz
        --species                    hsa or mmu (use hsa for human; mmu for mouse)
        --Trim_Mode                  Select a Small RNA-seq protocol
                1 NextFlex Small RNA-Seq Kit v3 (Hartwell Center protocl)
                2 NEBNext + sapkogrp/TrueSeq Index 26
                3 custom ; need to set adapter sequence by -A
        --help | h                    Show this usage statement.
       Note:
         All the above arguments can be configured by the command line interface or in the nextflow.config (default)


cat fq1.lst 
R022_3P /research/groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/R022_3P_L001_R1_001.fastq.gz
R044_3P /research/groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/R044_3P_L002_R1_001.fastq.gz

nextflow run sRNAtools.nf -profile local --fqlist fq1.lst --Trim_Mode 1 --species hsa --outdir sapkogrp --prefix sapkogrp 
```

## Output - example1

```bash
└── sapkogrp
    ├── clean
    │   ├── R022-12P
    │   │   └── R022-12P_clean.fa
    │   └── R022-3P
    │       └── R022-3P_clean.fa
    ├── results
    │   ├── R022-12P_table.result.txt
    │   ├── R022-3P_table.result.txt
    │   ├── sapkogrp_hsa_all.count.txt
    │   ├── sapkogrp_hsa_all.RPM.txt
    │   ├── sapkogrp_hsa_miRNA.count.txt
    │   └── sapkogrp_hsa_miRNA.RPM.txt
    └── Trim_Galore
        ├── R022-12P
        │   ├── R022-12P_R1.fastq.gz_trimming_report.txt
        │   ├── R022-12P_trimmed_fastqc.html
        │   └── R022-12P_trimmed.fq
        └── R022-3P
            ├── R022-3P_R1.fastq.gz_trimming_report.txt
            ├── R022-3P_trimmed_fastqc.html
            └── R022-3P_trimmed.fq
```

## sRNAtools
https://academic.oup.com/bib/article/22/1/463/5686255
http://rnainformatics.org.cn/sRNAtools/download.php

## Maintainers
* [Hongjian Jin]

[Hongjian Jin]: https://github.com/hongjianjin

## Update log
* 2023-04-18,  Version 1.0.0.  implement NEXTflex small RNAseq protocol
* 2023-04-20, Version 1.2.0. implement NEBNext small RNAseq protocol
* 2024-11-16, Version 1.3.0. check compatibility with rhel8; minor changes 
