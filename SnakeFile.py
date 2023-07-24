# A file to snakemake the workflow.
# Steps: Fastqc -> multiqc on fastqs 
# STAR on RNA-seq FASTQs
# make VCF with aidans script
# call variants - aidans
# call demuxlet - aidans
# scRNAseq pipeline
 



rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
