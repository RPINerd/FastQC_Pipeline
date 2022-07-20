# FastQC Pipeline
A straightforward program which can take a file list of fastq's, concatenate them, and perform QC on them in serial or scalable parallel

## Requirements
Python - built on 3.9x but should be fine for the forseeable future
FastQC - the underlying tool to perform QC on your read files
Bash Environment - current mechanics rely on system calls to the bash command line

## Input File
FastQC Pipe expects a tab-delimited text file with the following structure:

`Path/of/directory	|	Sample_Name	|	#_of_Lanes	|	R1/2/Both`

For example
```
/home/RPINerd/M01234/Fastq_Generation Exp001_S1   4   2
/home/RPINerd/M01234/Fastq_Generation Exp001_S2   4   1
/home/RPINerd/M01234/Fastq_Generation Exp001_S3   4   Both
```

## Running Pipe
General format for run is `python3 fastqc_pipe.py -f file_list.tsv`

The input file is the only required argument however there are some additional options
	- `-t #` / `--threads` specifies how many threads to allocate to the fastqc algorithm, currently capped at 12
	- `-m /home/RPINerd/merged/` / `--merge` allows for a preferred location to be specified for the merged lanes to be written to
	- `-v` / `--verbose` pumps out a ton of extra info and saves it to the file fastqc_pipe.log
