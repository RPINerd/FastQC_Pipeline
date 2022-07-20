import sys
import os
import subprocess
import logging
import argparse
from pathlib import Path

# FastQC pipeline       |       RPINerd, 07/20/22
# 
# FastQC_pipe.py will take an input of run files and sequentially analyze them
# with the fastqc tool. Input format is expected to be a tab-separated list 
# with the following columns:
# 
#   Path/of/directory   Sample_Name     #_of_Lanes  R1/2/Both
#
# ex.
#   /home/RPINerd/M01234/Fastq_Generation Exp001_S1   4   2
#   /home/RPINerd/M01234/Fastq_Generation Exp001_S2   4   1
#   /home/RPINerd/M01234/Fastq_Generation Exp001_S3   4   Both



# Sub to hunt down red oct.. I mean all the individual lane files for each readset
def collect_reads(rootpath, readset, lanes, r_int):

    matches = []
    read_match = readset + "_L00[1-" + lanes + "]_R" + r_int + "*.fastq*"

    logging.debug("read_match: " + read_match)

    for path in Path(rootpath).rglob(read_match):
        matches.append( str(path.resolve()).replace(" ", "\ ") )

    logging.debug("matches: " + str(matches))
    return matches


# Merge all the lanes individual files into a single fastq
def merge_fastq(r_int, r_files, sample_id):

    # Convert list into string for execution
    r_string = ' '.join(r_files)

    # Concatenate lanes into single file
    merge_name = sample_id + "_R" + r_int + ".fastq"
    #TODO check if compressed or not
    if r_string.find("gz") > 0:
        cmd = "zcat " + r_string + " > " + merge_name
    else:
        cmd = "cat " + r_string + " > " + merge_name

    logging.debug("merge_fastq\nr_string: " + r_string + "\nmerge_name: " + merge_name + "\ncmd: " + cmd + "\n")

    logging.info("Calling cat merge of Read " + str(r_int) + " files...")
    subprocess.run(cmd, shell=True)
    logging.info("Merge Complete!")
    return merge_name


# Just execute fastqc on given file
def process_files(file_list, threads):

    # Join array of files into string
    file_string = ' '.join(file_list)

    # Execute FastQC analysis
    #TODO handle both compressed and uncompressed
    fqc = "fastqc -t " + str(threads) + " " + file_string 
    subprocess.call(fqc, shell=True)
  

def main():

    # Argument Parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Your input *.tsv/*.csv with list of fastq files", required=True)
    parser.add_argument("-t", "--threads", help="Number of simultaneous threads to run", required=False, default=4, type=int)
    parser.add_argument("-m", "--merge", help="If desired, specify a location to save the fastq files after lane merge", required=False, default=False)
    parser.add_argument("-v", "--verbose", help="Outputs a lot more information for debugging and saves log", required=False, action='store_true')
    args = parser.parse_args()

    # Logging
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

    # Open up provided runlist file
    in_file = args.file
    if not os.path.isfile(in_file):
        logging.info('Error: input file does not exist!')
        return
    
    # Check for valid thread count
    # TODO how to get max number of threads available on machine?
    max_threads = 12
    if args.threads > max_threads:
        logging.info("Error: too many threads requested!")
        return
    else:
        threads = args.threads

    file_list = []
    # Parse input file, merge lanes, save reads into array for fastqc processing
    for run in open(in_file).read().splitlines():

        # Header line
        # TODO for future standardization, should either use bool flag for first run or designate first row character
        if run.startswith("#"):
            logging.info("Header Line...Skipping\n")
            continue

        cols = run.split('\t')

        path = cols[0]
        sample_id = cols[1]
        lanes = cols[2]
        reads = cols[3]

        logging.info("\n--Currently Processing--" \
            + "\nSample: " + sample_id \
            + "\nLanes: " + lanes \
            + "\nPath: " + path \
            + "\nReads: " + reads + "\n")

        if reads.upper() == "BOTH":
            for r in (1,2):

                r_file_list = collect_reads(path, sample_id, lanes, str(r))
                file_list.append(merge_fastq(str(r), r_file_list, sample_id))
        else:
            r_file_list = collect_reads(path, sample_id, lanes, str(reads))
            file_list.append(merge_fastq(str(reads), r_file_list, sample_id))

    # Pass the list of merged files to fastqc for processing
    process_files(file_list, threads)

    # Cleanup intermediates/logging
    if not args.verbose:
        os.remove("fastqc_pipe.log")
    if not args.merge:
        for file in file_list:
            os.remove(file)

if __name__ == "__main__":
    main()
    