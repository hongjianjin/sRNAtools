#!/usr/bin/env nextflow
/***********************************************************************************
 This is a Container of Nextflow pipeline-sRNAtools
 Input: Fastq filelist   
 Output: sRNA count tables
 Written by: Hongjian Jin, Ph.D.
 The Center for Applied Bioinformatics,St. Jude Children's Research Hosptial
 Date: 04/11/2023
 Last update: 04/18/2023
***********************************************************************************/
def DispConfig() {
log.info """
Welocme to run Nextflow Pipeline sRNAtools.nf [version 1.2.0, 4/20/2023]
Your configuration are the following:
  fqlist                : ${params.fastq_filelist}
  outdir                : ${params.outdir}
  prefix                : ${params.prefix}
  Trim_Mode	        : ${params.Trim_Mode}
  Select_Species        : ${params.Select_Species}
  A                     : ${params.A}
  The following detail information is based your configuration: 
 
"""  
//  exit 0    
}

def helpMessage() {
  log.info """
        Welocme to run Nextflow Pipeline sRNAtools.nf [version 1.2.0, 4/20/2023]
       Usage:
        A typical command for running the pipeline is as follows:
        nextflow run sRNAtools.nf -profile local --fqlist fq1.lst --Trim_Mode 1 --outdir run1 --prefix hendegrpq 
        nextflow run sRNAtools.nf -profile local --fqlist fq2.lst --Trim_Mode 2 --outdir run2 --prefix sapkogrp 

        Configurable arguments:
        --outdir                      The directory for the pipeline output
        --prefix                      output prefix for your analysis 
        --A			      adapter sequence for read 1 if Trim_Mode is 3  
        --fqlist                      Fastq file list that you want to call sRNAtools
		 format1 (2 columns) for single_end_oneLane:  
			  SampleName L001_R1_001.fastq.gz
		 format2 (2 columns) for single_end_multiLanes: 
			  SampleName L001_R1_001.fastq.gz,L002_R1_001.fastq.gz

	--Trim_Mode	              Select a Small RNA-seq protocol
                1 NextFlex Small RNA-Seq Kit v3 (Hartwell Center protocl)
                2 NEBNext + sapkogrp/TrueSeq Index 26
		3 custom ; need to set adapter sequence by -A
        --help | h                    Show this usage statement.
       Note:
         All the above arguments can be configured by the command line interface or in the nextflow.config (default)  
        """
}


// Show help message
if ((params.help) || (params.h)) {//|| (params.Trim_Mode) 
    helpMessage()
    exit 0
}


if ( params.Trim_Mode == 1){ // NEXTflex
    A="TGGAATTCTCGGGTGCCAAGG"
}else if ( params.Trim_Mode == 2){ // NEBNext? 
    A =" AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT -a TGGAATTCTCGGGTGCCAAGG -a GGAAGAGCACACGTCTGAACTCCAGTCACCAACTAATCTCGTATGCCGTC -n 4"
    //A =" AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT -a TGGAATTCTCGGGTGCCAAGG -a GGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTC -n 4"
    // R1 adapter:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    // R2 adapter: TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT to remove read through
    //GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT // RNA PCR Primer (RP1) 
    // TrueSeq index 26: GGAAGAGCACACGTCTGAACTCCAGTCACCAACTAATCTCGTATGCCGTC
    // trimming of IUPAC bases should work just fine
    //TrueSeq_Idx: GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCT
    //TrueSeq_Idx_trimed:     GGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTC
}else if ( params.Trim_Mode == 3){ // custom
     //Example for trimming multiple adapters
    //A =" AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TGGAATTCTCGGGTGCCAAGG -a GGAAGAGCACACGTCTGAACTCCAGTCACCAACTAATCTCGTATGCCGTC -n 3"
    A = params.A
    println ("custom adapter: $params.A  ")
}else{
    helpMessage()
    exit 0
}

// Display the configuration
DispConfig() 

/**********************************************************************************************
* This Process is responsible for trimming reads by Trim_Galore/0.6.6
***********************************************************************************************/

process Trim_Galore_NEXTflex {
   publishDir "${params.outdir}/Trim_Galore/${SampleName}", mode: 'copy', overwrite: true
   input:
   tuple val(SampleName), path(Read_Fastq) 
   
   output:
   tuple val(SampleName), path("${SampleName}_trimmed.fq") ,path("*.txt") ,path("*.html")  

   script:
   """
   trim_galore \
     --cores 4 \
     --output_dir . \
     --basename ${SampleName} \
     --dont_gzip \
     --fastqc_args "--outdir . " \
     --clip_R1 4 --three_prime_clip_R1 4 -a ${A} \
     ${Read_Fastq}
   """
}



/**********************************************************************************************
* This Process is responsible for trimming reads by Trim_Galore/0.6.6
***********************************************************************************************/
//https://github.com/FelixKrueger/TrimGalore/issues/86
process Trim_Galore_Custom {
   publishDir "${params.outdir}/Trim_Galore/${SampleName}", mode: 'copy', overwrite: true
   input:
   tuple val(SampleName), path(Read1_Fastq) 
   
   output:
   tuple val(SampleName), path("${SampleName}_trimmed.fq"),path("*.txt") ,path("*.html")  

   script:
   """
   trim_galore \
     --cores 4 \
     --output_dir . \
     --basename ${SampleName} \
     --dont_gzip \
     --fastqc_args "--outdir . " \
     -a "${A}" \
     ${Read1_Fastq}
   """
}


/*This process is responsible for zcat the multiple lanes' fastq.gz into one mergred fastq.gz  */
process Cat_SingleEndFastqGz_MultiLanes {
   publishDir "${params.outdir}/Fastq/${SampleName}", mode: 'symlink', overwrite: true
   input:
   tuple val(SampleName), path(Read1_Fastq)
   
   output:
   tuple val(SampleName), path("${SampleName}_R1.fastq.gz")

   script:
   """
     if [[ "$Read1_Fastq" == *","* ]]; then
	 echo "Info: cat multiple lanes."
	 cat `echo $Read1_Fastq | tr ',' ' '` >${SampleName}_R1.fastq.gz
     else
         ln -s $Read1_Fastq ${SampleName}_R1.fastq.gz
     fi
   """  
}


/***********************************************************************************************
*  This process is responsible for adapting the trim_galore output to fq2collapedFa.pl
************************************************************************************************/
process Collapse_FqToFa {
   publishDir "${params.outdir}/clean/${SampleName}", mode: 'copy', overwrite: true
   input:
   tuple val(SampleName), path(trimmed_fq) ,path("*.txt") ,path("*.html")  

   output:
   tuple val(SampleName), path("${SampleName}_clean.fa") 

   script:
   """
   perl /research/groups/cab/projects/automapper/common/hjin/bin/sRNAtools/fq2collapedFa.pl -i ${SampleName}_trimmed.fq -o ${SampleName}_clean.fa
   """
}
/***********************************************************************************************
*  This process is responsible for run_animal.pl step 
************************************************************************************************/

process Run_Animal {
   publishDir "${params.outdir}/results/", mode: 'copy', overwrite: true
   input:
   tuple val(SampleName), path(clean_fa) 
   output:
   tuple val(SampleName), path("${SampleName}_table.result.txt") 
   script:
   """
   perl /research/groups/cab/projects/automapper/common/hjin/bin/sRNAtools/run_animal.pl \
          -infile ${SampleName}_clean.fa \
	  -species ${params.Select_Species} \
	  -outdir . \
	  -jobid ${SampleName}
   """

}
/***********************************************************************************************
*  This process is responsible for creating symbolic links
************************************************************************************************/
process Copy_Results_Aadptor {
   input:
    tuple val(SampleName), path(table_result_txt) 
   output:
    path("${SampleName}")

   script:
   """
    ln -s $table_result_txt ${SampleName}
   """  
}
/***********************************************************************************************
*  This process is responsible for merging tables across samples
************************************************************************************************/
process Join_Tables {
   publishDir "${params.outdir}/results/", mode: 'copy', overwrite: true
   input:
   path (sRNA_result_list) 

   output:
   path("*.txt") 

   script:
   """
     Rscript  /research/groups/cab/projects/automapper/common/hjin/bin/sRNAtools/program/sRNAtools_merge_matrix.R \
     -L ${sRNA_result_list.collect{ "$it" }.join(",")} \
     -s ${params.Select_Species} \
     -O ${params.prefix}
   """
}

//****************************************************************************************************/

workflow {
    // Set up the one external input channel for samples (fastq files)  
    Read_Fastq_Ch  = Channel
	    .fromPath(params.fqlist)
	    .splitText()
	    .splitCsv(sep: '\t')

   if ( params.Trim_Mode == 1 ) { 
	Results_Ch=Cat_SingleEndFastqGz_MultiLanes(Read_Fastq_Ch)|Trim_Galore_NEXTflex|Collapse_FqToFa|Run_Animal|Copy_Results_Aadptor
   }else {
	Results_Ch=Cat_SingleEndFastqGz_MultiLanes(Read_Fastq_Ch)|Trim_Galore_Custom|Collapse_FqToFa|Run_Animal|Copy_Results_Aadptor
   }
    Join_Tables(Results_Ch.collect())

}

/*****************************************************************************************************
4/20/2023 Updates 
  1. delete params.project
  2. implement -Trim_mode 2/3
  3. use workflow and replace set with tuple
  4. update sRNAtools_merge_matrix.R to handle empty file ,single input absence of miRNA etc.
  5. implement Cat_SingleEndFastqGz_MultiLanes
/*****************************************************************************************************/