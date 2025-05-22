#!/usr/bin/env nextflow

/*
    Hisat2-Align aligns the input trimmed file and generates the bam output file.
    input data types: FASTQ, FASTA files
    output data type: BAM

    # reference fasta file download:
    https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz

    build index command:
    docker run --rm -v $PWD:/data/ gudmaprbk/hisat2.2.1:1.0.1 /usr/bin/hisat2-2.2.1/hisat2-build\
     /data/Homo_sapiens.GRCh38.dna_sm.toplevel.fa Homo_sapiens.GRCh38.dna_sm.toplevel.fa

    align command:
    docker run -v $PWD:/data/ gudmaprbk/hisat2.2.1:1.0.1 hisat2 -p <threads>\
     --summary-file <align_summary_file> -x <ref_file> -1 <read1_input_file>\
      -2 <read2_input_file> | samtools sort -@ 8 -o <output_file>

*/
params.align_outdir = "${PROJECT_DIR}/results/alignment"
params.ref_dir = "${PROJECT_DIR}/data"

process align_index {
    
    container 'gudmaprbk/hisat2.2.1:1.0.1'
    publishDir "${params.ref_dir}", mode: 'copy'
    label 'process_low'
    cpus 8

    input:
        path ref_fasta_file
    output:
        tuple val(ref_fasta_file.baseName), path("${ref_fasta_file.baseName}.*.ht2")
    script:
    """
    hisat2-build -p ${task.cpus} ${params.ref_dir}/${ref_file} ${ref_file}
    """
}


process align_single {

    tag {"Hisat2-align ${input_file.getBaseName(2)}"}
    label 'process_low'
    cpus 7

    container 'gudmaprbk/hisat2.2.1:1.0.1'
    publishDir "${params.align_outdir}", mode: 'copy'

    input:
        tuple path(input_file), path(ref_fasta_path)
    output:
        path("${input_file.getBaseName(2)}_alignsort.bam")
        // path "${sample_name}_align_summary.txt", emit: align_summary
    script:
    def sample_name = "${input_file.getBaseName(2)}"
    def align_outfile = "${sample_name}_alignsort.bam"
    """
    if [ -d ${params.align_outdir} ]; then
        continue
    else
        mkdir ${params.align_outdir}
    fi
    
    echo ${sample_name}
    hisat2 -p ${task.cpus} --summary-file ${sample_name}_align_summary.txt -x ${params.ref_dir}/${ref_fasta_path}\
     -U ${input_file} | samtools sort -@ ${task.cpus - 3} -o ${align_outfile}
    """
}


process align_paired {

    tag {"Hisat2-Align ${sample_name}"}
    label 'process_high'
    cpus 7

    container 'gudmaprbk/hisat2.2.1:1.0.1'
    publishDir "${params.align_outdir}", mode: 'copy'

    input:
        tuple val(sample_name), path(fq1_trimmed), path(fq2_trimmed), path(ref_fasta_path)
    output:
        path("${sample_name}_alignsort.bam")
    script:
    def align_outfile = "${sample_name}_alignsort.bam"
    """
    if [ -d ${params.align_outdir} ]; then
        continue
    else
        mkdir ${params.align_outdir}
    fi

    hisat2 -p ${task.cpus} --summary-file ${sample_name}_align_summary.txt -x ${params.ref_dir}/${ref_fasta_path}\
    -1 ${fq1_trimmed} -2 ${fq2_trimmed} | samtools sort -@ ${task.cpus - 3} -o ${align_outfile}
    """
}