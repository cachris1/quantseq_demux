"""
Utilities for running demuxlet.
"""
import logging
import os
import subprocess


def run_demuxlet_command(
    vcf_file: str,
    bam_file: str,
    barcode_file: str,
    output_directory: str,
    umi_tag: str,
    threads_per_job: int,
    threads: int = 1,
) -> None:

    """call demuxlet command.

    Args:
        vcf_file: vcf file
        bam_file: sorted bam files.
        barcode_file: directory to the barcode file.
        output_directory: output_directory
        umi_tag: umi tag UB or pN
        threads_per_job: number of threads per job
        threads: Number of threads to use during sorting.

    Returns:
        None
    """
    libName = bam_file.split("/")[-1].removesuffix(".srt.bam").removesuffix(".sort.bam")
    command_demuxlet = f"popscle demuxlet --sam {bam_file} --vcf {vcf_file}\
    --field PL \
    --out {output_directory}/demux_out_{libName} \
    --tag-UMI {umi_tag} \
    --group-list {barcode_file} \
    --min-callrate 0.2 \
    --alpha 0.0 --alpha 0.5"

    logging.info(f"calling demuxlet on files {bam_file} with {threads} threads.")
    os.system(command_demuxlet)