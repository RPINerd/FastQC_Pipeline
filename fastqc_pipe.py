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
import sys
import time
from pathlib import Path

#! Critical potential issue, if the order of lanes is ever inconsistent, everything breaks. Currently it
#! outputs a list in the order L004, L002, L003, L001, need to force order here to prevent headaches in
#! the future!


# --- Analysis --- #
# Sub to hunt down red oct.. I mean all the individual lane files for each readset
def collect_reads(rootpath, readset, readNumber):
    matches = []
    read_match = f"{readset}_L00[1-4]_R{readNumber}*.fastq*"

    for path in Path(rootpath).rglob(read_match):
        matches.append(str(path.resolve()).replace(" ", "\\ "))

    logging.debug(f"read_match regex:\t{read_match}\nmatches:\t{matches}")

    return matches


# Merge all the lanes individual files into a single fastq
def merge_fastq(jobs, merge_dir):
    merge_names = []
    processes = []
    for job in jobs:
        readNumber, readFiles, sample_id = job
        r_string = " ".join(readFiles)
        merge_name = f"{merge_dir}/{sample_id}_R{readNumber}.fastq"
        if r_string.find("gz"):
            merge_name += ".gz"
        merge_names.append(merge_name)

        # cmd = ["cat"]
        # cmd.extend(readFiles)
        # cmd.append(" > ")
        # cmd.append(merge_name)

        logging.debug(f"merge_fastq\nr_string:\t{r_string}\nmerge_name:\t{merge_name}\n")
        logging.info(f"Launching merge for {sample_id} R{readNumber}...")

        with open(merge_name, "wb") as concat:
            for file in readFiles:
                shutil.copyfileobj(open(file, "rb"), concat)
        # processes.append(subprocess.Popen(cmd, stdout=subprocess.PIPE))

    # still_running = True
    # total_procs = len(processes)
    # current_status = 0.00

    # while still_running:
    #     time.sleep(5)
    #     completed_procs = 0
    #     for proc in processes:
    #         logging.debug(proc.communicate())
    #         if proc.poll() is not None:
    #             completed_procs += 1
    #     new_status = round((completed_procs / total_procs) * 100, 2)
    #     if new_status != current_status:
    #         print(f"Merge status: {new_status}%", end="\r")
    #         current_status = new_status
    #     still_running = True if completed_procs < total_procs else False

    logging.info("Merge status: 100%\nMerge Completed!")
    logging.debug(f"Merge files final: {merge_names}")
    return merge_names


# Parse input file and collect all reads for each job
def parse_input_file(args):
    merge_jobs = []
    with open(args.file, "r") as runlist:
        logging.info("--Creating Merge Jobs--")
        for line in runlist:
            # Header line
            if line.startswith("#"):
                logging.info("Header Line...Skipping\n")
                continue

            sample_id = line.strip()
            reads = ["1", "2"] if args.reads == 3 else args.reads

            logging.info(f"Sample:\t{sample_id}\tReads:\t{reads}")

            for read in reads:
                read_file_list = collect_reads(args.dir, sample_id, read)
                if read_file_list == []:
                    logging.warning(f"No files were found for SampleID {sample_id}! Skipping...")
                else:
                    merge_jobs.append([read, read_file_list, sample_id])

    logging.info(f"{len(merge_jobs)} total jobs created.")

    return merge_jobs


# Just a tiny caller to fastqc
def fastqc_files(file_list, threads):
    fqc = ["fastqc", "-t", str(threads)].extend(file_list)
    subprocess.run(fqc, stdout=subprocess.PIPE)


# --- Utility --- #
def cli_parse():
    parser = argparse.ArgumentParser()
    input_type = parser.add_argument_group()
    # TODO allow inferring of reads by just providing a target folder
    input_type.add_argument("-d", "--dir", help="Directory where all fastq files are stored", required=True)
    input_type.add_argument("-f", "--file", help="Your input *.tsv/*.csv with list of fastq files", required=True)
    parser.add_argument(
        "-t", "--threads", help="Number of simultaneous threads to run", required=False, default=4, type=int
    )
    parser.add_argument(
        "-m",
        "--merge",
        help="If desired, specify a location to save the fastq files after lane merge",
        required=False,
        default="",
    )
    parser.add_argument(
        "-c", "--clean", help="After run, clean up the merge files from the disk", required=False, action="store_true"
    )
    parser.add_argument(
        "-r", "--reads", help="Choose whether to limit QC to only R1 or R2", required=False, choices=[1, 2], default=3
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
        logging.basicConfig(filename="fastqc_pipe.log", encoding="utf-8", level=getattr(logging, "DEBUG", None))
    else:
        logging.basicConfig(encoding="utf-8", level=getattr(logging, "INFO", None))

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(message)s")
    handler.setFormatter(formatter)
    root = logging.getLogger()
    root.addHandler(handler)


# --- Drivers --- #
def main(args) -> None:
    # Parse input for merge jobs
    merge_jobs = parse_input_file(args)

    # Merge all lanes into single file
    qc_jobs = merge_fastq(merge_jobs, args.merge)

    # Pass the list of merged files to fastqc for processing
    fastqc_files(qc_jobs, args.threads)

    # Cleanup intermediates/logging
    if args.clean:
        for file in qc_jobs:
            os.remove(file)


if __name__ == "__main__":
    # Parse user arguments and spin up logging
    args = cli_parse()
    setup_logging(args.verbose)
    logging.info("Logging started!")

    # Check for valid thread count
    max_threads = len(os.sched_getaffinity(0))
    logging.debug(f"Cores reported: {max_threads}")
    assert (
        args.threads <= max_threads
    ), f"Error: Too many threads requested! Maximum available on this machine is {max_threads}"

    # Validate runlist file
    assert os.path.isfile(args.file), f"Error: Input file ({args.file}) does not exist!"

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
