#!/usr/bin/env nextflow

// import modules
include { fastqc; fastqc_se } from "./modules/fastqc.nf"
include { trimmomatic_single; trimmomatic_paired } from "./modules/trimmomatic.nf"
include { align_index; align_single; align_paired } from "./modules/align.nf"
include {featurecounts} from "./modules/counts.nf"

// input data directory, input reference fasta and gtf files and adapter file path
params.data_dir = "${PROJECT_DIR}/data"
params.adapter_se_path = "${params.data_dir}/TruSeq2-SE.fa"
params.adapter_pe_path = "${params.data_dir}/TruSeq2-PE.fa"
params.ref_fasta_path = "${params.data_dir}/Homo_sapiens.GRCh38.dna_sm.toplevel.fa"
params.ref_gtf_file = "${params.data_dir}/Homo_sapiens.GRCh38.114.gtf"

// params.read_type = "PE"
params.read_type = "SE"
params.paired_end_files = "${params.data_dir}/reads/sample_*_{R1,R2}.f*.gz"
params.single_end_files = "${params.data_dir}/reads/sample_*_R1.f*.gz"


workflow {
    
    ref_fasta_ch = Channel.fromPath(params.ref_fasta_path, checkIfExists:true)
    ref_gtf_ch = Channel.fromPath(params.ref_gtf_file, checkIfExists:true)

    // Paired-End
    if (params.read_type == "PE") {

        // channel for paired input files
        Channel
            .fromFilePairs(params.paired_end_files)
            .ifEmpty {exit 1, "no fastq files found in the given path"}
            .view()
            .set{ read_pair_ch }

        // fastqc - pre trimming
        fastqc(read_pair_ch)
        
        // trimmomatic
        // adapter
        adapter_ch = Channel.fromPath(params.adapter_pe_path, checkIfExists: true)
        trimmomatic_paired(read_pair_ch.combine(adapter_ch))

        // Hisat2-align
        align_paired(trimmomatic_paired.out.combine(ref_fasta_ch))
        
        // featureCounts
        counts_combine_ch = align_paired.out.collect()
        featurecounts(counts_combine_ch, ref_gtf_ch)
    
    }
    else { // Single-End

        // Channel for input files (single-end)
        Channel
            .fromPath(params.single_end_files)
            .ifEmpty{exit 1, "no input fastq files found in the given path"}
            .view()
            .set{ read_ch }
        
        // fastqc - pre trimming
        fastqc_se(read_ch.collect())

        // trimmomatic
        adapter_ch = Channel.fromPath(params.adapter_se_path, checkIfExists: true)
        trimmomatic_single(read_ch.combine(adapter_ch))

        // Hisat2-align
        align_single(trimmomatic_single.out.combine(ref_fasta_ch))

        // featureCounts
        counts_combine_ch = align_single.out.collect()
        featurecounts(counts_combine_ch, ref_gtf_ch)
    }
}
