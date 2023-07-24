# A file to snakemake the workflow.
# Steps: Fastqc -> multiqc on fastqs 
# STAR on RNA-seq FASTQs
# make VCF with aidans script
# call variants - aidans
# call demuxlet - aidans
# scRNAseq pipeline
 

configfile: "config.yaml"

def get_fastqs(wildcards):
    return config["samples"][wildcards.sample]

    

rule all:
    input:
        "./prcessed_data/scRNA_counts_matrix"

rule fastqc:
    input:
        get_fastqs
    output:
        "./fastqc/"
    shell:
        "fastqc {input}"
        
rule multiqc:
    input:
       "./fastqc"
    output:
       "./multiqc"
    shell:
       "multiqc {input}"

       
rule star:
    input:
       config
