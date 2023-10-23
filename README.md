# quantseq_demux

A snakemake pipline for demuxing quantseq pools. Mainly a re-implementation of the protocol in the quantseq pool github, and uses their demux program

To demux quantseq pool data, use the 'bulk_to_bam' folder.  You must change the filepaths in the config.yaml to match your system, and run "snakemake --n-cores {n cores}"
You also must download the quantseq demux program: https://github.com/Lexogen-Tools/idemux as well as STAR, cutadapt, SAMTOOLS, and umi-tools

