"""
Utilities for processing BAM/SAM files.
"""
import logging
import os
import subprocess

import pysam


def sort_alignments(bam_file: str, output_directory: str, threads: int = 1) -> None:
    """Sorts alignments.

    Args:
        bam_file: BAM (or SAM) file input.
        output_directory: Output directory.
        threads: Number of threads to use during sorting.

    Returns:
        None
    """

    logging.info(f"Sorting {bam_file} with {threads} threads.")
    sorted_file_name = bam_file.replace(".bam", ".sorted.bam").split("/")[-1]

    pysam.sort(
        "-o",
        f"{output_directory}/{sorted_file_name}",
        "-@",
        str(threads),
        "-m",
        "7G",
        bam_file,
    )
    pysam.index(f"{output_directory}/{sorted_file_name}", "-@", str(threads))


def deduplicate_umis(bam_file: str, output_directory: str) -> None:
    """Deduplicate UMIs.

    Deduplicates UMIs with UMI-tools.

    Args:
        bam_file: Sorted bam file.
        output_directory: Output directory.

    Returns:
        None.
    """
    logging.info(f"Deduplicating {bam_file}")
    deduplicated_file_name = bam_file.replace(".sorted.bam", ".sorted.dedup.bam").split(
        "/"
    )[-1]

    process = subprocess.run(
        [
            "umi_tools",
            "dedup",
            "-I",
            bam_file,
            "-S",
            f"{output_directory}/{deduplicated_file_name}",
        ]
    )


def call_variants(
    bam_file: str,
    output_directory: str,
    reference_genome: str,
    depth_filter: int = 250,
) -> None:
    """Calls variants from a bam file.

    Uses bcftools to call variants on a bam file given a reference genome.

    Args:
        bam_file: Sorted, deduplicated bam file.
        output_directory: output directory.
        reference_genome: location of reference genome.
        depth_filter: Depth filter for calling variants.
        threads: Number of threads to use.

    Returns:
        None.
    """
    logging.info(f"Calling variants from {bam_file}")
    filtered_bcf = f'{output_directory}/{bam_file.replace(".sorted.dedup.bam", ".ind.vcf").split("/")[-1]}'

    command_pileup = (
        f"bcftools mpileup -d 100 -C 50 -Ou -f {reference_genome} {bam_file}"
    )
    command_call = f"bcftools call -Ou -mv"
    command_filter = f"bcftools filter -e '%QUAL<20 || DP>{depth_filter} || MAF<0.2'"

    process_variants = subprocess.call(
        command_pileup
        + " | "
        + command_call
        + " | "
        + command_filter
        + f" > {filtered_bcf}",
        shell=True,
    )

    process = subprocess.run(["bgzip", filtered_bcf])
    process = subprocess.run(["tabix", "-p", "vcf", f"{filtered_bcf}.gz"])
