import sys
import os
import subprocess
import logging
import argparse
import multiprocessing
from pathlib import Path

# FastQC pipeline       |       RPINerd, 01/26/23
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
def collect_reads(rootpath, readset, lanes, readNumber):

    matches = []
    #TODO can we just assume 4 lanes and just not append file if not found?
    #- Assert a file is found
    read_match = "{}_L00[1-{}]_R{}*.fastq*".format(readset, lanes, readNumber)

    for path in Path(rootpath).rglob(read_match):
        matches.append( str(path.resolve()).replace(" ", "\ ") )

    logging.debug("read_match:\t{}\nmatches:\t{}".format(read_match, matches))
    
    return matches


# Merge all the lanes individual files into a single fastq
def merge_fastq(readNumber, readFiles, sample_id):

    r_string = ' '.join(readFiles)

    merge_name = "{}_R{}.fastq".format(sample_id, readNumber)
    
    cat = "zcat" if r_string.find(".gz") else "cat"
    cmd = "{} {} > {}".format(cat, r_string, merge_name)

    logging.debug("merge_fastq\nr_string:\t{}\nmerge_name:\t{}\ncmd:\t{}\n".format(r_string, merge_name, cmd))

    logging.info("Merging R{} files...".format(readNumber))
    subprocess.run(cmd, shell=True)
    logging.info("Merge Complete!")
    
    return merge_name


def fastqc_files(file_list, threads):

    #TODO handle both compressed and uncompressed
    fqc = "fastqc -t {} {}".format(threads, ' '.join(file_list))
    subprocess.call(fqc, shell=True)
  

def main(args):

    # Validate runlist file
    assert os.path.isfile(args.file), 'Error: input file does not exist!'
    
    # Check for valid thread count
    max_threads = multiprocessing.cpu_count()
    threads = args.threads
    assert threads <= max_threads, "Error: too many threads requested!"

    file_list = []
    # Parse input file, merge lanes, save reads into array for fastqc processing
    with open(args.file, "r") as in_file:
        for line in in_file:

            # Header line
            if line.startswith("#"):
                logging.info("Header Line...Skipping\n")
                continue

            cols = line.split('\t')

            path = cols[0]
            sample_id = cols[1]
            lanes = cols[2]
            reads = cols[3]

            logging.info("\n--Currently Processing--\nSample:\t{}\nLanes:\t{}\nPath:\t{}\nReads:\t{}\n".format(sample_id, lanes, path, reads))

            #TODO Must be a cleaner way..
            if reads.upper() == "BOTH":
                for r in (1,2):

                    r_file_list = collect_reads(path, sample_id, lanes, str(r))
                    file_list.append(merge_fastq(str(r), r_file_list, sample_id))
            else:
                r_file_list = collect_reads(path, sample_id, lanes, str(reads))
                file_list.append(merge_fastq(str(reads), r_file_list, sample_id))

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
    #TODO allow inferring of reads by just providing a target folder
    parser.add_argument("-f", "--file", help="Your input *.tsv/*.csv with list of fastq files", required=True)
    parser.add_argument("-t", "--threads", help="Number of simultaneous threads to run", required=False, default=4, type=int)
    parser.add_argument("-m", "--merge", help="If desired, specify a location to save the fastq files after lane merge", required=False, default=False)
    parser.add_argument("-v", "--verbose", help="Outputs a lot more information for debugging and saves log", required=False, action='store_true')
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
    