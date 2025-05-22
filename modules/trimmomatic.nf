#!/usr/bin/env nextflow

/*
    Trimmomatic module for trimming the input data.
    input data types: FASTQ, FASTA files

    command:
    time docker run -v $PWD:/data/ staphb/trimmomatic trimmomatic <read_type:SE/PE> -threads <threads>\
     -trimlog <trim_log_file> -summary <trim_summary_file> <input_files> <output_files>\
      ILLUMINACLIP:<adapter_file>:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
*/
params.trim_outdir = "${PROJECT_DIR}/results/trimmomatic"

process trimmomatic_single {

    tag {"Trimmomatic ${input_file}"}
    label 'process_low'
    cpus 8

    container 'staphb/trimmomatic:0.39'
    publishDir "${params.trim_outdir}", mode: 'copy'

    input:
        tuple path(input_file), path(adapter_file)
    
    output:
        path("${input_file.getBaseName(2)}_trimmed.fastq.gz")
        // path "${input_file.getBaseName(2)}_trim_summary.txt", emit: trim_summary
    
    script:
    def sample_name = "${input_file.getBaseName(2)}"
    """
    if [ -d ${params.trim_outdir} ]; then
        continue
    else
        mkdir ${params.trim_outdir}
    fi

    echo ${sample_name}
    trimmomatic SE -threads ${task.cpus} -trimlog ${sample_name}_trimlog.txt -summary ${sample_name}_trim_summary.txt ${input_file}\
    ${sample_name}_trimmed.fastq.gz ILLUMINACLIP:${adapter_file}:2:30:10\
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process trimmomatic_paired {

    tag {"Trimmomatic ${sample_name}"}
    label 'process_low'
    queue 'normal'
    cpus 8

    container 'staphb/trimmomatic:0.39'
    publishDir "${params.trim_outdir}", mode: 'copy'

    input:
        tuple val(sample_name), path(input_file), path(adapter_file)
    output:
        tuple val(sample_name), path("${sample_name}_R1_trimmed.fastq.gz"), path("${sample_name}_R2_trimmed.fastq.gz")
        // path "${sample_name}_trim_summary.txt", emit: trim_summary
    script:
    def (fq1, fq2) = input_file
    def fq1_trimmed_read = "${sample_name}_R1_trimmed.fastq.gz"
    def fq1_untrimmed_read = "${sample_name}_R1_untrimmed.fastq.gz"
    def fq2_trimmed_read = "${sample_name}_R2_trimmed.fastq.gz"
    def fq2_untrimmed_read = "${sample_name}_R2_untrimmed.fastq.gz"
    """
    if [ -d ${params.trim_outdir} ]; then
        continue
    else
        mkdir ${params.trim_outdir}
    fi

    trimmomatic PE -threads ${task.cpus} -trimlog ${sample_name}_trimlog.txt -summary ${sample_name}_trim_summary.txt ${fq1} ${fq2}\
    ${fq1_trimmed_read} ${fq1_untrimmed_read} ${fq2_trimmed_read} ${fq2_untrimmed_read} ILLUMINACLIP:${adapter_file}:2:30:10\
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
