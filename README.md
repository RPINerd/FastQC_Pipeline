# FastQC Pipeline

A straightforward program which can take a file list of fastq's, concatenate them, and perform FastQC on them in serial or scalable parallel

## Requirements

[Python](https://www.python.org/) - built on 3.10x but should be fine for the forseeable future  
[FastQC](https://github.com/s-andrews/FastQC) - the underlying tool to perform QC on your read files  
Bash Environment - current mechanics rely on system calls to the bash command line  

## Input File

File structure should follow the convention `SampleID_Lane#_Read#.fastq`. For example:
Exp001_S1_L002_R1_001.fastq

FastQC Pipe expects a tab-delimited text file with the following structure (lines led with '#' are ignored):

`Path/of/directory | Sample_Name | R1/2/Both`

Example:

```csv
/home/RPINerd/M01234/Fastq_Generation Exp001_S1   2
/home/RPINerd/M01234/Fastq_Generation Exp001_S2   1
/home/RPINerd/M01234/Fastq_Generation Exp001_S3   Both
```

## Running Pipeline

General format for run is `python3 fastqc_pipe.py -f file_list.tsv`

The input file is the only required argument however there are some additional options:  
`--threads` or `-t #` specifies how many threads to allocate to the fastqc algorithm, currently capped at 12  
`--merge` or `-m <desired/path/to/>` allows for a preferred location to be specified for the merged lanes to be written to  
`--verbose` or `-v` pumps out a ton of extra info and saves it to the file fastqc_pipe.log  
