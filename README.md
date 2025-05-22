# RNASeq_Pipeline_nextflow
RNASeq Pipeline, from raw reads to counts matrix by using open source bioinformatic tools such as FastQC, Trimmomatic, Hisat2-align and featureCounts with Nextflow, workflow language.

To running this RNASeq pipeline,
   - To keep all the input files in data/reads directory.
   - keep all the reference files such as fasta file and hisat2-index files, gtf file in /data directory.
   - keep the adapter file also in /data directory.

To mention the read type paired-end or single-end like, in the rnaseq.nf file
```
params.read_type = "SE"
(or)
params.read_type = "PE"
```
else, mention the read type info in the command line itself.
```
nextflow run rnaseq.nf -profile docker --read_type SE
```

To mention the input files in the rnaseq.nf file like
```
params.paired_end_files = "${params.data_dir}/reads/sample_*_{R1,R2}.f*.gz" for paired-end
params.single_end_files = "${params.data_dir}/reads/sample_*_R1.f*.gz" for single-end
```

To run the pipeline for single-end input files:
```
docker run -it --rm -v /var/run/docker.sock:/var/run/docker.sock -v /usr/bin/docker:/usr/bin/docker -e PROJECT_DIR=$PWD\
 -v $PWD:$PWD -w $PWD nextflow/nextflow:24.10.5 nextflow run rnaseq.nf -profile docker --read_type SE
```
for Paired-end input files:
```
docker run -it --rm -v /var/run/docker.sock:/var/run/docker.sock -v /usr/bin/docker:/usr/bin/docker -e PROJECT_DIR=$PWD\
 -v $PWD:$PWD -w $PWD nextflow/nextflow:24.10.5 nextflow run rnaseq.nf -profile docker --read_type PE
```
