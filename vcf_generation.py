"""
Main logic behind reference VCF generation.

This file stores the main entry point for the VCF generation pipeline, and makes
heavy use of the high level functionality in other modules of the
vevo-pipeline codebase.
"""
import argparse
import configparser
import logging
import os
import subprocess
from functools import partial

import ngs_tools as ngs
from joblib import delayed

from vevotx.mixins import logger
from vevotx.utilities import (bam_utilities, logging_utilities,
                              pipeline_setup_utilities)


@logging_utilities.log_all_args
def make_reference_vcf(
    bam_directory: str,
    output_directory: str,
    reference_genome: str,
    threads: int,
    min_threads_per_job: int,
    depth_filter: int,
    verbose: bool,
) -> None:

    # itemize the sample BAMs to be processed.
    sample_bams = []
    for _file in os.listdir(bam_directory):
        if _file.endswith(".bam") or _file.endswith(".sam"):
            sample_bams.append(_file)

    logger.info(f"Processing the following BAMs: {sample_bams}")

    # compute number of jobs that can be handled by compute environment
    n_jobs_parallel = (
        len(sample_bams)
        if (min_threads_per_job * len(sample_bams)) <= threads
        else (threads // min_threads_per_job)
    )
    threads_per_job = threads // n_jobs_parallel
    print(threads_per_job, n_jobs_parallel)

    # 1. sort each bam --------------------------------------------------------
    sort_partial = partial(bam_utilities.sort_alignments, threads=threads_per_job)
    ngs.utils.ParallelWithProgress(
        n_jobs=n_jobs_parallel,
        total=len(sample_bams),
        desc="Sorting alignments",
    )(
        delayed(sort_partial)(f"{bam_directory}/{sample_bam}", output_directory)
        for sample_bam in sample_bams
    )

    # 2. Deduplicate UMIs -----------------------------------------------------
    deduplicate_partial = partial(bam_utilities.deduplicate_umis)
    sorted_bam_files = [
        bam_file
        for bam_file in os.listdir(output_directory)
        if bam_file.endswith("sorted.bam")
    ]
    ngs.utils.ParallelWithProgress(
        n_jobs=threads,
        total=len(sample_bams),
        desc="Deduplicating UMIs",
    )(
        delayed(deduplicate_partial)(
            f"{output_directory}/{sorted_bam_file}", output_directory
        )
        for sorted_bam_file in sorted_bam_files
    )

    # 3. Call Variants -----------------------------------------------------
    call_variants_partial = partial(
        bam_utilities.call_variants,
        depth_filter=depth_filter,
    )
    deduplicated_bam_files = [
        bam_file
        for bam_file in os.listdir(output_directory)
        if bam_file.endswith("sorted.dedup.bam")
    ]
    ngs.utils.ParallelWithProgress(
        n_jobs=threads,
        total=len(deduplicated_bam_files),
        desc="Calling variants",
    )(
        delayed(call_variants_partial)(
            f"{output_directory}/{deduplicated_bam_file}",
            output_directory,
            reference_genome,
        )
        for deduplicated_bam_file in deduplicated_bam_files
    )

    # 4. Merge VCFs -----------------------------------------------------
    print([vcf for vcf in os.listdir(output_directory) if vcf.endswith("ind.vcf.gz")])
    with open(f"{output_directory}/reference.vcf", "w") as reference_out:
        subprocess.run(
            ["vcf-merge"]
            + [
                f"{output_directory}/{vcf}"
                for vcf in os.listdir(output_directory)
                if vcf.endswith("ind.vcf.gz")
            ],
            stdout=reference_out,
        )

    logging.info(
        f"Done with pipeline. Reference written to {output_directory}/reference.vcf"
    )


@logger.namespaced("vcf-generation-main")
@logging_utilities.log_runtime
def main():

    # --------------- Create Argument Parser & Read in Arguments -------------- #
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Specify a config file for analysis.")

    args = parser.parse_args()

    config_filepath = args.config

    minimum_parameters = [
        "output_directory",
        "reference_genome",
        "threads",
        "bam_directory",
    ]

    with open(config_filepath, "r") as f:
        pipeline_parameters = pipeline_setup_utilities.parse_config(
            f.read(), minimum_parameters
        )

    # pull out parameters
    output_directory = pipeline_parameters["vcf-generation"]["output_directory"]
    bam_directory = pipeline_parameters["vcf-generation"]["bam_directory"]
    reference_genome = pipeline_parameters["vcf-generation"]["reference_genome"]
    min_threads_per_job = pipeline_parameters["vcf-generation"]["min_threads_per_job"]
    depth_filter = pipeline_parameters["vcf-generation"]["depth_filter"]
    threads = pipeline_parameters["vcf-generation"]["threads"]
    logger_name = pipeline_parameters["vcf-generation"]["log_path"]
    verbose = pipeline_parameters["vcf-generation"].get("verbose", False)
    if verbose:
        logger.setLevel(logging.DEBUG)

    # set up output directory
    pipeline_setup_utilities.setup(
        output_directory, logger_name=logger_name, verbose=verbose
    )

    # start up pipeline
    make_reference_vcf(
        bam_directory,
        output_directory,
        reference_genome,
        threads,
        min_threads_per_job,
        depth_filter,
        verbose,
    )


if __name__ == "__main__":
    main()
