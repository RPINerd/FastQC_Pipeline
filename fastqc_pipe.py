"""
    FastQC Pipeline | RPINerd, 12/09/24

    FastQC_pipe.py will take an input of run files and analyze them with the fastqc tool in a quasi-parallel mode.
    Input format is expected to be a list with the just the read file IDs:

    Exp001_S1
    Exp001_S2
    Exp001_S3
"""

import argparse
import logging
import os
import shutil
import subprocess
from pathlib import Path


def arg_parser() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()

    # TODO allow inferring of reads by just providing a target folder
    input_type = parser.add_argument_group()
    input_type.add_argument(
        "-d",
        "--dir",
        help="Directory where all fastq files are stored",
        required=True,
    )
    input_type.add_argument(
        "-f",
        "--file",
        help="Your input *.tsv/*.csv with list of fastq files",
        required=True,
    )

    parser.add_argument(
        "-t",
        "--threads",
        help="Number of simultaneous threads to run. By default it will use 1 thread per sample, \
            or all available threads. Whichever is lower.",
        required=False,
        default=0,
        type=int,
    )
    parser.add_argument(
        "-m",
        "--merge",
        help="If desired, specify a location to save the fastq files after lane merge",
        required=False,
        default="",
    )
    parser.add_argument(
        "-c",
        "--clean",
        help="After run, clean up the merge files from the disk",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--reads",
        help="Choose whether to limit QC to only R1 or R2. By defaut QC is run on both.",
        required=False,
        choices=[1, 2],
        default=3,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Outputs a lot more information for debugging and saves log",
        required=False,
        action="store_true",
    )

    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    """
    Setup logging for the script

    Args:
        verbose (bool): If True, log all the things

    Returns:
        None
    """
    if verbose:
        logging.basicConfig(
            filename="fastqc_pipe.log",
            filemode="a",
            format="%(asctime)s - %(levelname)s - %(message)s",
            encoding="utf-8",
            datefmt="%H:%M:%S",
            level=logging.DEBUG,
        )
    else:
        logging.basicConfig(
            format="%(asctime)s - %(message)s",
            encoding="utf-8",
            datefmt="%M:%S",
            level=logging.INFO,
        )


def collect_reads(rootpath: str, readset: str, read_number: str) -> list:
    """
    Sub to hunt down red oct.. I mean all the individual lane files for each readset

    Args:
        rootpath (str): Root directory to search for files
        readset (str): Sample ID to search for
        read_number (str): Read number to search for

    Returns:
        matches (list): List of all files found matching the readset and read_number
    """
    matches = []
    read_match = f"{readset}_L00[1-4]_R{read_number}*.fastq*"

    for path in Path(rootpath).rglob(read_match):
        matches.append(str(path.resolve()).replace(" ", "\\ "))
    matches.sort()

    logging.debug(f"read_match regex:\t{read_match}\nmatches:\t{matches}")
    return matches


def merge_fastq(jobs: list[tuple], merge_dir: str) -> list:
    """
    Merge all the lanes individual files into a single fastq

    Args:
        jobs (list): List of jobs to merge
        merge_dir (str): Directory to save the merged files

    Returns:
        merge_names (list): List of all the merged files
    """
    logging.info("Beginning Lane Files Merge...")
    merge_names = []
    for job in jobs:
        read_number, read_files, sample_id = job
        r_string = " ".join(read_files)
        merge_name = f"{sample_id}_R{read_number}.fastq"
        if merge_dir:
            merge_name = merge_dir + "/" + merge_name
        if r_string.find("gz"):
            merge_name += ".gz"
        merge_names.append(merge_name)

        logging.info(f"Merging {sample_id} R{read_number}...")

        # TODO must test and handle non-zipped fastq files
        with Path.open(merge_name, "wb") as concat:
            for file in read_files:
                shutil.copyfileobj(Path.open(file, "rb"), concat)
                logging.info(f"Merge: {str(file).split('/')[-1]} -> {merge_name}")
        logging.info(f"Done {merge_name}")

    logging.info("Lane Files Merge Completed!")
    logging.debug(f"Merge files final: {merge_names}")
    return merge_names


def parse_input_file(args: argparse.Namespace) -> list:
    """
    Parse input file and collect all reads for each job

    Args:
        args (argparse.Namespace): Parsed arguments from the user

    Returns:
        merge_jobs (list): List of all jobs to be merged
    """
    merge_jobs = []
    with Path.open(args.file) as runlist:
        logging.info("Parsing sample list...")
        for line in runlist:
            # Header line
            if line.startswith("#"):
                logging.info("Header Line...Skipping\n")
                continue

            sample_id = line.strip()
            reads = ["1", "2"] if args.reads == 3 else args.reads

            logging.debug(f"Sample:\t{sample_id}\tReads:\t{reads}")

            for read in reads:
                read_file_list = collect_reads(args.dir, sample_id, read)
                if read_file_list == []:
                    logging.warning(f"No files were found for SampleID {sample_id}! Skipping...")
                else:
                    merge_jobs.append([read, read_file_list, sample_id])

    logging.info(f"{len(merge_jobs)} total jobs created.")

    return merge_jobs


def main(args: argparse.Namespace) -> None:
    """Main function to run the FastQC pipeline"""
    # Parse input for merge jobs
    merge_jobs = parse_input_file(args)

    # Merge all lanes into single file
    qc_jobs = merge_fastq(merge_jobs, args.merge)

    # Establish number of threads to use for FastQC
    threads = args.threads
    if threads == 0:
        # args.threads=0 means auto-detect and either match threads to jobs or max out available threads
        threads = len(qc_jobs) if len(qc_jobs) <= len(os.sched_getaffinity(0)) else len(os.sched_getaffinity(0))

    # Pass the list of merged files to fastqc for processing
    fqc = ["fastqc", "-t", str(threads)]
    fqc.extend(qc_jobs)
    logging.debug(f"FastQC Command: {fqc}")
    subprocess.run(fqc, stdout=subprocess.PIPE, check=False)

    # Cleanup intermediates/logging
    if args.clean:
        for file in qc_jobs:
            Path.unlink(file)


if __name__ == "__main__":
    # Parse user arguments and spin up logging
    args = arg_parser()
    setup_logging(args.verbose)
    logging.info("Logging started!")

    # Validate runlist file
    assert Path.is_file(args.file), f"Error: Input file ({args.file}) does not exist!"

    # Check for valid thread count
    max_threads = len(os.sched_getaffinity(0))
    logging.debug(f"Cores reported: {max_threads}")
    if args.threads > max_threads:
        logging.warning(f"Too many threads requested! Maximum available on this machine is {max_threads}.")

    # Create the desired merge directory if needed
    if args.merge:
        Path(args.merge).mkdir(parents=True, exist_ok=True)

    # Check for FastQC install
    app = shutil.which("fastqc")
    logging.debug(f"Shutil reports app as {app}")
    if app is not None:
        try:
            o = subprocess.check_output([app, "-h"], stderr=subprocess.STDOUT)
        except ChildProcessError:
            raise "FastQC application was not found!"

    # Execute Pipeline
    main(args)
