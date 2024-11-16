## CAB small RNA-seq Nextflow pipeline [V1.2.0] based on sRNAtools

## CAB internal wikipage
https://wiki.stjude.org/display/CABI/CAB+miRNAseq+nextflow+pipeline

## Dependencies 
* conda environment:  /research/groups/cab/projects/automapper/common/hjin/bin/miniconda3_envs/sRNAtools 
* databases: /research/groups/cab/projects/automapper/common/tchang1/Softwares/sRNAtools/db/  (defined by DBCONFIG.txt) 
* bowtie-1.2.2 index: /research/groups/cab/projects/automapper/common/hjin/bin/bowtie-1.2.2/indexes
* samtools/1.9 , bedtools/2.25.0, bowtie/1.2.2
* trimgalore/0.6.6,fastqc/0.11.8, python/3.7.0
* R/3.6.1
* nextflow/21.10.5

## Input

* format1 (2 columns) for single_end_oneLane: 
  SampleName R1.fastq.gz 
* format2 (2 columns) for single_end_multiLanes: 
  SampleName L1_R1.fastq.gz,L2_R1.fastq.gz 

## Run - example1
```bash
conda activate /research/groups/cab/projects/automapper/common/hjin/bin/miniconda3_envs/sRNAtools
export PATH="/research/groups/cab/projects/automapper/common/hjin/bin/miniconda3_envs/sRNAtools/bin:${PATH}"
export R_LIBS=/research/rgs01/applications/hpcf/authorized_apps/rhel7_apps/R/install/3.6.1/lib64/R/library/
nextflow run sRNAtools.nf --help 

cat fq1.lst 
R022_3P /research/groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/R022_3P_L001_R1_001.fastq.gz
R044_3P /research/groups/sapkogrp/projects/CAB/common/UMD_sRNA_seq/R044_3P_L002_R1_001.fastq.gz

nextflow run sRNAtools.nf -profile local --fqlist fq1.lst --Trim_Mode 1 --outdir sapkogrp --prefix sapkogrp 
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

