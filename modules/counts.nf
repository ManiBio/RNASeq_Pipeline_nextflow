#!/usr/bin/env nextflow

/*
    FeatureCounts generates the gene count matrix file from alignment BAM file.
    input:
        BAM
    output:
        featurecounts.csv, featurecounts.csv.summary

    # Reference GTF file download
    wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz

    featurecounts command with default arguments:
    docker run --rm -v $PWD:/data/ -w /data/ pegi3s/feature-counts:latest featureCounts -T 8 -t exon -g gene_id\
     -a <ref_gtf_file> -o <featurecounts_counts.csv file>  <input_bam files>

*/

params.counts_outdir = "${PROJECT_DIR}/results/counts"


process featurecounts {

    tag {"featureCounts"}
    label 'process_low'
    cpus 4

    container 'pegi3s/feature-counts:2.0.0'
    publishDir "${params.counts_outdir}", mode: 'copy'

    input:
        path(alignedsort_reads)
        path(ref_gtf_path)

    output:
        path "featureCounts.csv", emit: counts_file
        path "featureCounts_filtered.csv", emit: filtered_counts_file
    script:
    """
    if [ -d ${params.counts_outdir} ]; then
        continue
    else
        mkdir ${params.counts_outdir}
    fi

    featureCounts -T ${task.cpus} -t exon -g gene_id -a ${ref_gtf_path} -o featureCounts.csv ${alignedsort_reads}
    grep -v '#' featureCounts.csv | cut -f1,7- > featureCounts_filtered.csv
    """
}