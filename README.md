# Long Read Structural Variations for Oysters
A collection of executed scripts for Crassostrea gigas PacBio reads.


## Brief Background
A suite of scripts utilised for processing the Pacific oyster PacBio reads is available.


## Citation
Hyungtaek Jung et al.: Unraveling the biotechnological potential of Crassostrea gigas: comparative genomics & structural variations, [BioRxiv](https://www.biorxiv.org/XXXX).


## Contents:
+ [LICENSE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#license) 
+ [GETTING STARTED](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#getting-started)
+ [FAQ](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#faq)
+ [WIKI PAGE](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#wiki-page)
+ [AUTHORS](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#authors)
+ [COPYRIGHT](https://github.com/OZTaekOppa/FastQHandler/blob/main/README.md#copyright)



## LICENSE
The current scripts are available under the MIT license and incorporate various open-source software. 


## GETTING STARTED
The present collection of scripts is primarily designed for analysing PacBio reads of the Pacific oyster.


### Tested Files
- C. gigas1: Genome assembly [GCA_011032805.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_011032805.1), Project number [PRJNA598006](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA598006) and Raw data [SRA_598006](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=598006) in NCBI
- C. gigas2: Genome assembly [GCF_902806645.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_902806645.1), Project number [PRJEB35351](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB35351) and Raw data [SRA_596972](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=596972) in NCBI


### General Usage
Each program has provided its own description.


### fqread_stats (qrs)
- To generate fastq statistics for a fastq file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A summary of fastq with its diverse statistics (e.g. sequence length, max length, N50, and more).

Example usage
```
python3 fqread_stats.py --input-seq test_dna.fq --out output_stats.csv --t 1 --mem 2
```
+ It will accept both short and long-read FASTQ files.

- Parameter explanation
	1. python 3: Call python 3
	1. fqread_stats.py:  Call fqread_stats module
	1. python3 fqread_stats.py --help: Check help menu
		+ --input-seq: Indicate an input multi-line fastq file and its path
		+ --out: Indicate an output csv file with the summary of statistics
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fqread_stats_unlim (qsu)
- To generate fastq statistics for multiple fastq files
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Multiple fastq files.
	+ Output: A summary of fastq with its diverse statistics (e.g. sequence length, max length, N50, and more).

Example usage
```
python3 fqread_stats_unlim.py --input-seqs test_dna1.fq test_dna2.fq test_dna3.fq test_dna4.fq test_dna5.fq --out ./ --t 1 --mem 2
```
+ It will accept both short and long-read FASTQ files.

- Parameter explanation
	1. python 3: Call python 3
	1. fqread_stats_unlim.py:  Call fqread_stats_unlim module
	1. python3 fqread_stats_unlim.py --help: Check help menu
		+ --input-seqs: Indicate multiple fastq files and theirs paths
		+ --out: Indicate an output folder to save a csv file with the summary of statistics
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### fq2fa_convert (qac)
- To convert a fastq file as a fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fastq file.
	+ Output: A converted fasta file with the same fastq header.

Example usage
```
python3 fq2fa_convert.py --input-seq test_dna.fq --out output_convert.fq --t 1 --mem 2
```
+ It will accept both short and long-read FASTQ files.

- Parameter explanation
	1. python 3: Call python 3
	1. fq2fa_convert.py:  Call fq2fa_convert module
	1. python3 fq2fa_convert.py --help: Check help menu
		+ --input-seq: Indicate an input fastq file and its path
		+ --out: Indicate a converted output fasta file with the same fastq header
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)





## FAQ

We encourage users to use the [Issues](https://github.com/OZTaekOppa/lrsv_oysters/issues).


## WIKI PAGE

Please see the GitHub page.


## AUTHORS

**Hyungtaek Jung** and et al 2024.


## COPYRIGHT

The full scripts are distributed under the MIT license. 
