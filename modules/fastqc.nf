#!/usr/bin/env nextflow

/*
    FastQC module for preprocessing the input data.
    input data types: FASTQ, FASTA, SAM and BAM
*/
params.fastqc_outdir = "${PROJECT_DIR}/results/fastqc"

process fastqc {

    tag {"FASTQC ${sample_name}"}
    label 'process_low' 
    queue = 'normal'
    cpus 2

    container 'staphb/fastqc:latest'
    publishDir "${params.fastqc_outdir}", mode: 'copy'

    input:
        tuple val(sample_name), path(input_file)
    
    output:
        tuple val(sample_name), path("*.html"), path("*.zip")
    script:
    """
    echo ${sample_name}
    if [ -d ${params.fastqc_outdir} ]; then
        continue
    else
        mkdir ${params.fastqc_outdir}
    fi
    fastqc -t ${task.cpus} ${input_file}
    """
}

process fastqc_se {

    label 'process_low' 
    queue = 'normal'
    cpus 2

    container 'staphb/fastqc:latest'
    publishDir "${params.fastqc_outdir}", mode: 'copy'

    input:
        path(input_file)
    
    output:
        tuple path("*.html"), path("*.zip")
    script:
    """
    if [ -d ${params.fastqc_outdir} ]; then
        continue
    else
        mkdir ${params.fastqc_outdir}
    fi

    fastqc -t ${task.cpus} ${input_file}
    """
}