import argparse
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

'''
    FastQC Pipeline | RPINerd, 03/08/23

    FastQC_pipe.py will take an input of run files and analyze them with the fastqc tool in a quasi-parallel mode.
    Input format is expected to be a tab-separated list with the following columns:

    Path/of/directory   Sample_Name     Read1/2/Both

    /home/RPINerd/M01234/Fastq_Generation Exp001_S1   2
    /home/RPINerd/M01234/Fastq_Generation Exp001_S2   1
    /home/RPINerd/M01234/Fastq_Generation Exp001_S3   Both
  '''

# Sub to hunt down red oct.. I mean all the individual lane files for each readset
def collect_reads(rootpath, readset, readNumber):

    matches = []
    read_match = f"{readset}_L00[1-4]_R{readNumber}*.fastq*"

    for path in Path(rootpath).rglob(read_match):
        matches.append( str(path.resolve()).replace(" ", "\ ") )

    logging.debug(f"read_match regex:\t{read_match}\nmatches:\t{matches}")

    return matches if len(matches) else -1


# Merge all the lanes individual files into a single fastq
def merge_fastq(jobs):

    merge_names = []
    processes = []
    for job in jobs:

        readNumber, readFiles, sample_id = job
        r_string = ' '.join(readFiles)
        merge_name = f"{sample_id}_R{readNumber}.fastq"
        merge_names.append(merge_name)

        cat = "cat" if r_string.find("gz") == -1 else "zcat"
        cmd = f"{cat} {r_string} > {merge_name}"

        logging.debug(f"merge_fastq\nr_string:\t{r_string}\nmerge_name:\t{merge_name}\ncmd:\t{cmd}\n")
        logging.info(f"Launching merge for {sample_id} R{readNumber}...")

        processes.append(subprocess.Popen(cmd, shell=True))

    still_running = True
    total_procs = len(processes)

    while still_running:

        time.sleep(2)
        status = 0
        for proc in processes:
            if proc.poll() is not None:
                status += 1

        print(f"Merge status: {round((status/total_procs)*100,2)}%", end="\r")
        still_running = True if status < total_procs else False

    logging.info("Merge status: 100%\nMerge Completed!")
    return merge_names


# Parse input file and collect all reads for each job
def parse_input_file(file):

    merge_jobs = []
    with open(file, "r") as runlist:

        logging.info("--Creating Merge Jobs--")
        for line in runlist:

             # Header line
            if line.startswith("#"):
                logging.info("Header Line...Skipping\n")
                continue

            cols = line.strip().split('\t')

            path = cols[0]
            sample_id = cols[1]
            reads = str(cols[2])
            reads = ['1', '2'] if reads.lower() == "both" else reads

            logging.info(f"Sample:\t{sample_id}\tReads:\t{reads}")

            for read in reads:
                read_file_list = collect_reads(path, sample_id, read)
                if read_file_list == -1:
                        logging.warning(f"No files were found for SampleID {sample_id}! Skipping...")
                        continue
                merge_jobs.append([read, read_file_list, sample_id])

    logging.info(f"{len(merge_jobs)} total jobs created.")

    return merge_jobs


# Just a tiny caller to fastqc
def fastqc_files(file_list, threads):

    #TODO handle both compressed and uncompressed
    files = ' '.join(file_list)
    fqc = f"fastqc -t {threads} {files}"
    subprocess.call(fqc, shell=True)


def main(args):

    # Parse input for merge jobs
    merge_jobs = parse_input_file(args.file)

    # Merge all lanes into single file
    qc_jobs = merge_fastq(merge_jobs)

    # Pass the list of merged files to fastqc for processing
    fastqc_files(qc_jobs, args.threads)

    # Cleanup intermediates/logging
    if not args.verbose:
        os.remove("fastqc_pipe.log")
    if args.clean:
        for file in qc_jobs:
            os.remove(file)


if __name__ == "__main__":

    # Argument Parsing
    parser = argparse.ArgumentParser()
    input_type = parser.add_argument_group()
    # TODO allow inferring of reads by just providing a target folder
    input_type.add_argument("-d", "--dir", help="Directory where all fastq files are stored", required=True)
    input_type.add_argument("-f", "--file", help="Your input *.tsv/*.csv with list of fastq files", required=True)
    parser.add_argument("-t", "--threads", help="Number of simultaneous threads to run", required=False, default=4, type=int)
    parser.add_argument("-m", "--merge", help="If desired, specify a location to save the fastq files after lane merge", required=False, default=False)
    parser.add_argument("-c", "--clean", help="After run, clean up the merge files from the disk", required=False, action='store_true')
    # TODO allow single input for R1/R2/Both that will apply to all reads
    #parser.add_argument("-r", "--reads", help="Choose whether to QC R1, R2 or both", required=False, choices=[1,2,""])
    parser.add_argument("-v", "--verbose", help="Outputs a lot more information for debugging and saves log", required=False, action='store_true')
    args = parser.parse_args()

    # Logging Setup
    if args.verbose:
        logging.basicConfig(filename='fastqc_pipe.log', encoding='utf-8', level=getattr(logging, "DEBUG", None))
    else:
        logging.basicConfig(encoding='utf-8', level=getattr(logging, "INFO", None))
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(message)s")
    handler.setFormatter(formatter)
    root = logging.getLogger()
    root.addHandler(handler)
    logging.info("Logging started!")

    # Check for valid thread count
    max_threads = len(os.sched_getaffinity(0))
    logging.debug(f"Cores reported: {max_threads}")
    assert args.threads <= max_threads, f"Error: Too many threads requested! Maximum available on this machine is {max_threads}"

    # Validate runlist file
    assert os.path.isfile(args.file), f"Error: Input file ({args.file}) does not exist!"

    # Check for FastQC install
    app = shutil.which("fastqc")
    logging.debug(f"Shutil reports app as {app}")
    if app != None:
        try:
            o = subprocess.check_output([app, '-h'],stderr= subprocess.STDOUT)
        except:
            raise "FastQC application was not found!"

    # Execute Pipeline
    main(args)

