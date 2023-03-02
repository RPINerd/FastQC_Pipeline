import argparse
import asyncio
import logging
import multiprocessing
import os
from pathlib import Path
import subprocess
import sys

# FastQC pipeline       |       RPINerd, 01/27/23
# 
# FastQC_pipe.py will take an input of run files and sequentially analyze them
# with the fastqc tool. Input format is expected to be a tab-separated list 
# with the following columns:
# 
#   Path/of/directory   Sample_Name     Read1/2/Both
#
# ex.
#   /home/RPINerd/M01234/Fastq_Generation Exp001_S1   2
#   /home/RPINerd/M01234/Fastq_Generation Exp001_S2   1
#   /home/RPINerd/M01234/Fastq_Generation Exp001_S3   Both


# Sub to hunt down red oct.. I mean all the individual lane files for each readset
def collect_reads(rootpath, readset, readNumber):

    matches = []
    read_match = f"{readset}_L00[1-4]_R{readNumber}*.fastq*"

    for path in Path(rootpath).rglob(read_match):
        matches.append( str(path.resolve()).replace(" ", "\ ") )

    logging.debug(f"read_match:\t{read_match}\nmatches:\t{matches}")
    
    return matches if len(matches) else -1


# Merge all the lanes individual files into a single fastq
def merge_fastq(readNumber, readFiles, sample_id):

    r_string = ' '.join(readFiles)

    merge_name = f"{sample_id}_R{readNumber}.fastq"
    
    cat = "cat" if r_string.find("gz") == -1 else "zcat"
    cmd = f"{cat} {r_string} > {merge_name}"

    logging.debug(f"merge_fastq\nr_string:\t{r_string}\nmerge_name:\t{merge_name}\ncmd:\t{cmd}\n")

    logging.info(f"Merging R{readNumber} lanes...")
    subprocess.run(cmd, shell=True)
    logging.info("Merge Complete!")
    
    return merge_name


def fastqc_files(file_list, threads):

    #TODO handle both compressed and uncompressed
    files = ' '.join(file_list)
    fqc = f"fastqc -t {threads} {files}"
    subprocess.call(fqc, shell=True)
  

def main(args):

    # Validate runlist file
    assert os.path.isfile(args.file), 'Error: input file does not exist!'
    
    #TODO Check for fastqc install?
    
    # Check for valid thread count
    max_threads = multiprocessing.cpu_count()
    threads = args.threads
    assert threads <= max_threads, f"Error: too many threads requested! Maximum available on this machine is {max_threads}"

    file_list = []
    # Parse input file, merge lanes, save reads into array for fastqc processing
    with open(args.file, "r") as in_file:
        for line in in_file:

            # Header line
            if line.startswith("#"):
                logging.info("Header Line...Skipping\n")
                continue

            cols = line.strip().split('\t')

            path = cols[0]
            sample_id = cols[1]
            reads = str(cols[2])

            logging.info(f"\n--Currently Processing--\nSample:\t{sample_id}\nPath:\t{path}\nReads:\t{reads}\n")
            
            #TODO Must be a cleaner way..
            if reads.lower() == "both":
                for read in ['1', '2']:
                    r_file_list = collect_reads(path, sample_id, read)
                    if r_file_list == -1:
                        logging.warning(f"No files were found for SampleID {sample_id}! Skipping.")
                        continue
                    file_list.append(merge_fastq(read, r_file_list, sample_id))

            else:
                r_file_list = collect_reads(path, sample_id, reads)
                if r_file_list == -1:
                    logging.warning(f"No files were found for SampleID {sample_id}! Skipping.")
                    continue
                file_list.append(merge_fastq(reads, r_file_list, sample_id))

    # Pass the list of merged files to fastqc for processing
    fastqc_files(file_list, threads)

    # Cleanup intermediates/logging
    if not args.verbose:
        os.remove("fastqc_pipe.log")
    if not args.merge:
        for file in file_list:
            os.remove(file)


if __name__ == "__main__":

    # Argument Parsing
    parser = argparse.ArgumentParser()
    # TODO allow inferring of reads by just providing a target folder
    # TODO allow single input for R1/R2/Both that will apply to all reads
    parser.add_argument("-f", "--file", help="Your input *.tsv/*.csv with list of fastq files", required=True)
    parser.add_argument("-t", "--threads", help="Number of simultaneous threads to run", required=False, default=4, type=int)
    parser.add_argument("-m", "--merge", help="If desired, specify a location to save the fastq files after lane merge", required=False, default=False)
    parser.add_argument("-v", "--verbose", help="Outputs a lot more information for debugging and saves log", required=False, action='store_true')
    parser.add_argument("-h", "--help", help="Display the help dialog", required=False)
    args = parser.parse_args()

    # Logging Setup
    if args.verbose:
        logging.basicConfig(filename='fastqc_pipe.log', encoding='utf-8', level=getattr(logging, "DEBUG", None))
    else:
        logging.basicConfig(filename='fastqc_pipe.log', encoding='utf-8', level=getattr(logging, "INFO", None))
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    root = logging.getLogger()
    root.addHandler(handler)
    logging.info('Logging started!')

    # Execute Pipeline

    main(args)
    
