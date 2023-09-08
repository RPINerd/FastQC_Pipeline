"""
    FastQC Pipeline | RPINerd, 09/08/23

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


def cli_parse() -> argparse.ArgumentParser:
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
    args = parser.parse_args()

    return args


def setup_logging(verbose) -> None:
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


# Sub to hunt down red oct.. I mean all the individual lane files for each readset
def collect_reads(rootpath, readset, readNumber):
    matches = []
    read_match = f"{readset}_L00[1-4]_R{readNumber}*.fastq*"

    for path in Path(rootpath).rglob(read_match):
        matches.append(str(path.resolve()).replace(" ", "\\ "))
    matches.sort()

    logging.debug(f"read_match regex:\t{read_match}\nmatches:\t{matches}")
    return matches


# Merge all the lanes individual files into a single fastq
def merge_fastq(jobs, merge_dir):
    logging.info("Beginning Lane Files Merge...")
    merge_names = []
    for job in jobs:
        readNumber, readFiles, sample_id = job
        r_string = " ".join(readFiles)
        merge_name = f"{sample_id}_R{readNumber}.fastq"
        if merge_dir != "":
            merge_name = merge_dir + "/" + merge_name
        if r_string.find("gz"):
            merge_name += ".gz"
        merge_names.append(merge_name)

        logging.info(f"Merging {sample_id} R{readNumber}...")

        # TODO must test and handle non-zipped fastq files
        with open(merge_name, "wb") as concat:
            for file in readFiles:
                shutil.copyfileobj(open(file, "rb"), concat)
                logging.info(f"Merge: {str(file).split('/')[-1]} -> {merge_name}")
        logging.info(f"Done {merge_name}")

    logging.info("Lane Files Merge Completed!")
    logging.debug(f"Merge files final: {merge_names}")
    return merge_names


# Parse input file and collect all reads for each job
def parse_input_file(args):
    merge_jobs = []
    with open(args.file, "r") as runlist:
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


def main(args) -> None:
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
    subprocess.run(fqc, stdout=subprocess.PIPE)

    # Cleanup intermediates/logging
    if args.clean:
        for file in qc_jobs:
            os.remove(file)


if __name__ == "__main__":
    # Parse user arguments and spin up logging
    args = cli_parse()
    setup_logging(args.verbose)
    logging.info("Logging started!")

    # Validate runlist file
    assert os.path.isfile(args.file), f"Error: Input file ({args.file}) does not exist!"

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
